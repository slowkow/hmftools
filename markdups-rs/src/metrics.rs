//! picard.sam.DuplicationMetrics and the metrics file as htsjdk's MetricsFile writes it.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use hmf_common::fmt::{Rounding, decimal_format_rounding, java_double_to_string};

use crate::bam::Rec;

#[derive(Clone, Debug, Default)]
pub struct DuplicationMetrics {
    pub library: Vec<u8>,
    pub unpaired_reads_examined: i64,
    pub read_pairs_examined: i64,
    pub secondary_or_supplementary_rds: i64,
    pub unmapped_reads: i64,
    pub unpaired_read_duplicates: i64,
    pub read_pair_duplicates: i64,
    pub read_pair_optical_duplicates: i64,
    pub percent_duplication: Option<f64>,
    pub estimated_library_size: Option<i64>,
}

impl DuplicationMetrics {
    pub fn new(library: &[u8]) -> Self {
        DuplicationMetrics { library: library.to_vec(), ..Default::default() }
    }

    /// addReadToLibraryMetrics
    pub fn add_read(&mut self, rec: &Rec) {
        if rec.is_unmapped() {
            self.unmapped_reads += 1;
        } else if rec.is_secondary_or_supplementary() {
            self.secondary_or_supplementary_rds += 1;
        } else if !rec.is_paired() || rec.mate_unmapped() {
            self.unpaired_reads_examined += 1;
        } else {
            self.read_pairs_examined += 1;
        }
    }

    /// addDuplicateReadToMetrics
    pub fn add_duplicate(&mut self, rec: &Rec) {
        if !rec.is_secondary_or_supplementary() && !rec.is_unmapped() {
            if !rec.is_paired() || rec.mate_unmapped() {
                self.unpaired_read_duplicates += 1;
            } else {
                self.read_pair_duplicates += 1;
            }
        }
    }

    /// calculateDerivedFields
    pub fn calculate_derived_fields(&mut self) {
        self.estimated_library_size =
            estimate_library_size(self.read_pairs_examined - self.read_pair_optical_duplicates, self.read_pairs_examined - self.read_pair_duplicates);
        let examined = self.unpaired_reads_examined + self.read_pairs_examined;
        self.percent_duplication = Some(if examined != 0 {
            (self.unpaired_read_duplicates + self.read_pair_duplicates * 2) as f64 / (self.unpaired_reads_examined + self.read_pairs_examined * 2) as f64
        } else {
            0.0
        });
    }

    /// calculateRoiHistogram: None when there is no library size
    pub fn roi_histogram(&self) -> Option<Vec<(f64, f64)>> {
        let size = self.estimated_library_size?;
        let unique = self.read_pairs_examined - self.read_pair_duplicates;
        let mut h = Vec::with_capacity(100);
        let mut x = 1.0f64;
        while x <= 100.0 {
            h.push((x, estimate_roi(size, x, self.read_pairs_examined, unique)));
            x += 1.0;
        }
        Some(h)
    }
}

fn f(x: f64, c: f64, n: f64) -> f64 {
    c / x - 1.0 + (-n / x).exp()
}

/// DuplicationMetrics.estimateLibrarySize: the Lander-Waterman solution by bisection
pub fn estimate_library_size(read_pairs: i64, unique_read_pairs: i64) -> Option<i64> {
    let dups = read_pairs - unique_read_pairs;
    if read_pairs > 0 && dups > 0 {
        let mut m = 1.0f64;
        let mut big_m = 100.0f64;
        let (u, n) = (unique_read_pairs as f64, read_pairs as f64);
        // the IllegalStateException branch cannot be reached: f(c, c, n) = exp(-n / c) > 0
        while f(big_m * u, u, n) > 0.0 {
            big_m *= 10.0;
        }
        for _ in 0..40 {
            let r = (m + big_m) / 2.0;
            let v = f(r * u, u, n);
            if v == 0.0 {
                break;
            } else if v > 0.0 {
                m = r;
            } else if v < 0.0 {
                big_m = r;
            }
        }
        Some((u * (m + big_m) / 2.0) as i64)
    } else {
        None
    }
}

/// DuplicationMetrics.estimateRoi
pub fn estimate_roi(estimated_library_size: i64, x: f64, pairs: i64, unique_pairs: i64) -> f64 {
    let s = estimated_library_size as f64;
    s * (1.0 - (-(x * pairs as f64) / s).exp()) / unique_pairs as f64
}

/// FormatUtil.format(double): at most 6 fraction digits, HALF_DOWN, no grouping, "?" for NaN and
/// infinity.
pub fn format_double(v: f64) -> String {
    if v.is_nan() {
        return "?".into();
    }
    if v.is_infinite() {
        return if v < 0.0 { "-?".into() } else { "?".into() };
    }
    decimal_format_rounding(v, 6, Rounding::HalfDown)
}

/// The set-size histograms of LibraryIdGenerator: Histogram<Double> keyed by set size.
#[derive(Default, Clone, Debug)]
pub struct SetSizeHistograms {
    pub all_sets: BTreeMap<u64, f64>,
    pub optical_sets: BTreeMap<u64, f64>,
    pub non_optical_sets: BTreeMap<u64, f64>,
}

impl SetSizeHistograms {
    pub fn merge(&mut self, o: &SetSizeHistograms) {
        for (a, b) in [(&mut self.all_sets, &o.all_sets), (&mut self.optical_sets, &o.optical_sets), (&mut self.non_optical_sets, &o.non_optical_sets)] {
            for (k, v) in b {
                *a.entry(*k).or_insert(0.0) += v;
            }
        }
    }
}

/// MetricsFile.write: the two string headers, the metrics rows, then the histograms side by side.
pub fn write_metrics(command_line: &str, started: &str, metrics: &[DuplicationMetrics], roi: Option<&[(f64, f64)]>, sets: &SetSizeHistograms) -> String {
    let mut out = String::new();
    out.push_str("## htsjdk.samtools.metrics.StringHeader\n");
    let _ = writeln!(out, "# {command_line}");
    out.push_str("## htsjdk.samtools.metrics.StringHeader\n");
    let _ = writeln!(out, "# Started on: {started}");
    out.push('\n');
    if !metrics.is_empty() {
        out.push_str("## METRICS CLASS\tpicard.sam.DuplicationMetrics\n");
        out.push_str("LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n");
        for m in metrics {
            let _ = writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                String::from_utf8_lossy(&m.library),
                m.unpaired_reads_examined,
                m.read_pairs_examined,
                m.secondary_or_supplementary_rds,
                m.unmapped_reads,
                m.unpaired_read_duplicates,
                m.read_pair_duplicates,
                m.read_pair_optical_duplicates,
                m.percent_duplication.map(format_double).unwrap_or_default(),
                m.estimated_library_size.map(|v| v.to_string()).unwrap_or_default()
            );
        }
    }
    out.push('\n');
    // the histograms that have bins, in the order added: CoverageMult (one library only), all_sets,
    // optical_sets, non_optical_sets
    let mut columns: Vec<(&str, BTreeMap<OrdF64, f64>)> = Vec::new();
    let mut bin_label = "set_size";
    if let Some(r) = roi {
        if !r.is_empty() {
            columns.push(("CoverageMult", r.iter().map(|&(k, v)| (OrdF64(k), v)).collect()));
            bin_label = "BIN";
        }
    }
    for (label, h) in [("all_sets", &sets.all_sets), ("optical_sets", &sets.optical_sets), ("non_optical_sets", &sets.non_optical_sets)] {
        if !h.is_empty() {
            columns.push((label, h.iter().map(|(&k, &v)| (OrdF64(k as f64), v)).collect()));
        }
    }
    if !columns.is_empty() {
        let mut keys: Vec<OrdF64> = columns.iter().flat_map(|(_, h)| h.keys().copied()).collect();
        keys.sort();
        keys.dedup();
        out.push_str("## HISTOGRAM\tjava.lang.Double\n");
        out.push_str(bin_label);
        for (label, _) in &columns {
            out.push('\t');
            out.push_str(label);
        }
        out.push('\n');
        for k in keys {
            out.push_str(&java_double_to_string(k.0));
            for (_, h) in &columns {
                out.push('\t');
                out.push_str(&format_double(h.get(&k).copied().unwrap_or(0.0)));
            }
            out.push('\n');
        }
    }
    out.push('\n');
    out
}

#[derive(Clone, Copy, PartialEq, Debug)]
struct OrdF64(f64);
impl Eq for OrdF64 {}
impl PartialOrd for OrdF64 {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for OrdF64 {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&o.0)
    }
}

/// java.util.Date.toString() of the current time ("Sat Nov 28 04:19:37 GMT 2026"), in the local time
/// zone (GMT where none is set, as a JVM in a container without one).
pub fn java_date_now() -> String {
    let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs() as i64).unwrap_or(0);
    let mut tm: libc_tm = unsafe { std::mem::zeroed() };
    let t = now as libc_time_t;
    unsafe { localtime_r(&t, &mut tm) };
    let zone = if tm.tm_zone.is_null() { String::new() } else { unsafe { std::ffi::CStr::from_ptr(tm.tm_zone) }.to_string_lossy().to_string() };
    let zone = if zone.is_empty() || zone == "UTC" { "GMT".to_string() } else { zone };
    const DAYS: [&str; 7] = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"];
    const MONTHS: [&str; 12] = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
    format!(
        "{} {} {:02} {:02}:{:02}:{:02} {} {}",
        DAYS[tm.tm_wday.clamp(0, 6) as usize],
        MONTHS[tm.tm_mon.clamp(0, 11) as usize],
        tm.tm_mday,
        tm.tm_hour,
        tm.tm_min,
        tm.tm_sec,
        zone,
        tm.tm_year + 1900
    )
}

#[allow(non_camel_case_types)]
type libc_time_t = i64;

#[repr(C)]
#[allow(non_camel_case_types)]
struct libc_tm {
    tm_sec: i32,
    tm_min: i32,
    tm_hour: i32,
    tm_mday: i32,
    tm_mon: i32,
    tm_year: i32,
    tm_wday: i32,
    tm_yday: i32,
    tm_isdst: i32,
    tm_gmtoff: i64,
    tm_zone: *const std::ffi::c_char,
}

unsafe extern "C" {
    fn localtime_r(t: *const libc_time_t, tm: *mut libc_tm) -> *mut libc_tm;
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Picard's DuplicationMetricsTest / MarkDuplicates outputs: the library size and ROI of the small
    /// synthetic case run through GATK 4.6.1.0 (3 pairs, 2 duplicates, 1 optical)
    #[test]
    fn library_size_and_roi_as_gatk() {
        let mut m = DuplicationMetrics::new(b"lib1");
        m.unpaired_reads_examined = 2;
        m.read_pairs_examined = 3;
        m.unpaired_read_duplicates = 1;
        m.read_pair_duplicates = 2;
        m.read_pair_optical_duplicates = 1;
        m.calculate_derived_fields();
        assert_eq!(m.estimated_library_size, Some(1));
        assert_eq!(format_double(m.percent_duplication.unwrap()), "0.625");
        let roi = m.roi_histogram().unwrap();
        let got: Vec<String> = roi.iter().take(5).map(|&(_, v)| format_double(v)).collect();
        assert_eq!(got, ["0.950213", "0.997521", "0.999877", "0.999994", "1"]);
    }

    #[test]
    fn library_size_none_without_duplicates() {
        assert_eq!(estimate_library_size(100, 100), None);
        assert_eq!(estimate_library_size(0, 0), None);
        // DuplicationMetrics.main READ_PAIRS READ_PAIR_DUPLICATES, GATK 4.6.1.0's jar
        assert_eq!(estimate_library_size(10000, 9000), Some(46607));
        assert_eq!(estimate_library_size(1000000, 350000), Some(376419));
        assert_eq!(estimate_library_size(74000000, 24000000), Some(25373351));
    }

    #[test]
    fn half_down_formatting() {
        assert_eq!(format_double(0.0), "0");
        assert_eq!(format_double(1.0), "1");
        assert_eq!(format_double(0.6875), "0.6875");
        assert_eq!(format_double(123456789.0), "123456789");
        assert_eq!(format_double(f64::NAN), "?");
        assert_eq!(format_double(1.0 / 3.0), "0.333333");
        assert_eq!(format_double(2.0 / 3.0), "0.666667");
    }
}
