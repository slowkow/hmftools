//! VCF record text: genotype fields as VariantContextFactory.createGenotype + htsjdk encoding.

use crate::constants::*;
use crate::evidence::ReadContextCounter;
use crate::fmt::*;

pub const FORMAT_KEYS: &str =
    "GT:AD:AED:AF:AFBQ:AFMQ:ARCBQ:ARMQ:ASBQ:DP:FSB:MUC:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RSB:SAC:UMI_CNT";

const CORE_AF_MIN: i32 = 3;
const CORE_AF_FULL_RATIO: f64 = 0.25;

fn ints(v: &[i32]) -> String {
    v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}

pub fn genotype(c: &ReadContextCounter) -> String {
    let q = &c.qual_counters;
    let mut depth = c.depth();
    let mut alt_support = c.alt_support();
    let mut ref_support = c.counts.ref_;
    let strong = c.strong_alt_support();
    let avg = |total: f64, n: i32| -> i32 {
        if n > 0 {
            jround(total / n as f64) as i32
        } else {
            0
        }
    };
    let avg_map_quality = avg(q.map_quality_total as f64, depth);
    let avg_alt_map_quality = avg(q.alt_map_quality_total as f64, alt_support);
    let avg_base_quality = avg(q.recalibrated_base_quality_total as f64, depth);
    let avg_alt_base_quality = jround(c.average_alt_recalibrated_base_quality()) as i32;
    let avg_alt_final_base_quality = avg(q.alt_final_base_quality_total, strong);
    let avg_alt_final_map_quality = avg(q.alt_final_map_quality_total, strong);

    let mut vaf = c.vaf();
    let core_proportion = c.counts.core as f64 / c.counts.alt_support() as f64;
    if c.counts.core >= CORE_AF_MIN && core_proportion >= CORE_AF_FULL_RATIO {
        alt_support -= c.counts.core;
        vaf = alt_support as f64 / c.counts.total as f64;
    }
    if c.uncertain_core_base_count > 0 {
        let alt_ratio = alt_support as f64 / depth as f64;
        let ref_ratio = ref_support as f64 / depth as f64;
        depth += c.uncertain_core_base_count;
        ref_support += jround(c.uncertain_core_base_count as f64 * ref_ratio) as i32;
        alt_support += jround(c.uncertain_core_base_count as f64 * alt_ratio) as i32;
    }
    let umi = match c.consensus_type_counts {
        Some(v) => v,
        None => [c.depth(), 0, 0, c.counts.strong_support(), 0, 0],
    };
    let fields = [
        "./.".to_string(),
        ints(&[ref_support, alt_support]),
        format!(
            "{},{}",
            vcf_double(c.read_edge_distance.avg_distance_from_edge()),
            vcf_double(c.read_edge_distance.avg_alt_distance_from_edge())
        ),
        vcf_double(vaf),
        avg_alt_final_base_quality.to_string(),
        avg_alt_final_map_quality.to_string(),
        ints(&[avg_base_quality, avg_alt_base_quality]),
        ints(&[avg_map_quality, avg_alt_map_quality]),
        (c.average_alt_seq_tech_base_quality() as i32).to_string(),
        depth.to_string(),
        format!(
            "{},{}",
            java_fixed(c.non_alt_fragment_strand_bias.bias(), 3),
            java_fixed(c.alt_fragment_strand_bias.bias(), 3)
        ),
        c.fragment_coords.min_count().to_string(),
        ints(&c.counts.to_array()),
        c.improper_pair_count.to_string(),
        ints(&[c.jitter.shortened, c.jitter.lengthened]),
        ints(&c.qualities.to_array()),
        format!(
            "{},{}",
            java_fixed(c.non_alt_read_strand_bias.bias(), 3),
            java_fixed(c.alt_read_strand_bias.bias(), 3)
        ),
        c.simple_alt_matches.to_string(),
        ints(&umi),
    ];
    fields.join(":")
}

use crate::filters::{SageVariant, ALL_FILTERS};
use anyhow::{Context, Result};
use rust_htslib::bgzf;
use std::io::Write as IoWrite;

const HEADER_FIELDS: &str = include_str!("../resources/vcf_header_fields.txt");
pub const SAGE_VERSION: &str = "5.1";

/// One VCF data line, as VariantContextFactory.create + htsjdk's encoder would write it.
pub fn record(chr: &str, v: &SageVariant) -> String {
    let c = &v.candidate;
    let ctx = &c.ctx;
    let t = &v.tumor_counters[0];
    let mut info: Vec<(&str, Option<String>)> = Vec::new();
    if v.tier() == crate::types::Tier::Hotspot {
        info.push(("HOTSPOT", None));
    } else if v.near_hotspot {
        info.push(("NEAR_HOTSPOT", None));
    }
    if v.has_local_phase_sets() {
        info.push(("LPS", Some(ints(&t.local_phase_sets))));
        info.push(("LPS_RC", Some(ints(&t.lps_counts))));
    }
    info.push(("MED", Some(vcf_double(t.read_edge_distance.max_alt_distance_from_edge))));
    if let Some(h) = &ctx.homology {
        let bases = String::from_utf8_lossy(&h.bases).into_owned();
        info.push(("MH", Some(bases.clone())));
        info.push(("RC_MH", Some(bases)));
    }
    info.push(("MQF", Some(vcf_double(t.map_qual_factor))));
    if v.mixed_impact > 0 {
        info.push(("MSG", Some(v.mixed_impact.to_string())));
    }
    if v.near_indel {
        info.push(("NEARBY_INDEL", None));
    }
    info.push((
        "RC_INFO",
        Some(format!(
            "{}-{}-{}-{}-{}-{}",
            ctx.alignment_start,
            ctx.var_index,
            ctx.left_flank_str(),
            ctx.core_str(),
            ctx.right_flank_str(),
            ctx.read_cigar
        )),
    ));
    info.push(("RC_NM", Some(c.min_events.to_string())));
    if let Some(r) = &ctx.max_repeat {
        info.push(("RC_REPC", Some(r.count.to_string())));
        info.push(("RC_REPS", Some(r.bases_str().to_string())));
    }
    if let Some(r) = &ctx.max_ref_repeat {
        info.push(("REP_C", Some(r.count.to_string())));
        info.push(("REP_S", Some(r.bases_str().to_string())));
    }
    info.push(("TIER", Some(v.tier().name().to_string())));
    info.push(("TNC", Some(String::from_utf8_lossy(ctx.trinucleotide()).into_owned())));
    info.push(("TQP", Some(vcf_double(t.tumor_qual_probability))));
    info.sort_by(|a, b| a.0.cmp(b.0));
    let info_str = info
        .iter()
        .map(|(k, val)| match val {
            Some(x) => format!("{k}={x}"),
            None => k.to_string(),
        })
        .collect::<Vec<_>>()
        .join(";");

    // QUAL: log10PError = round(logTqp * 10) / 10, written as -10 x that with trailing ".00" removed
    let tqp = jmin(jmax(t.tumor_qual_probability, 0.0), 1.0);
    let log_tqp = jmax(tqp, 1e-20).log10();
    let log10_p_error = jround(log_tqp * 10.0) as f64 / 10.0;
    let qual = log10_p_error * -10.0 + 0.0;
    let mut qual_str = java_fixed(qual, 2);
    if let Some(stripped) = qual_str.strip_suffix(".00") {
        qual_str = stripped.to_string();
    }

    let filter_str = if v.filters.is_empty() {
        "PASS".to_string()
    } else {
        let mut names: Vec<&str> = v.filters.iter().map(|f| f.name()).collect();
        names.sort();
        names.join(";")
    };

    let mut line = format!(
        "{}\t{}\t.\t{}\t{}\t{}\t{}\t{}\t{}",
        chr,
        v.position(),
        v.variant().ref_str(),
        v.variant().alt_str(),
        qual_str,
        filter_str,
        info_str,
        FORMAT_KEYS
    );
    for rc in v.ref_counters.iter().chain(v.tumor_counters.iter()) {
        line.push('\t');
        line.push_str(&genotype(rc));
    }
    line
}

pub struct VcfWriter {
    path: String,
    out: bgzf::Writer,
    /// the previous region's records, held back as VcfWriter.checkLaterVariants may add to them
    pending: Vec<(i32, String)>,
}

impl VcfWriter {
    pub fn create(path: &str, contigs: &[(String, i32)], reference_ids: &[String], tumor_ids: &[String]) -> Result<Self> {
        let mut out = bgzf::Writer::from_path(path).with_context(|| format!("creating {path}"))?;
        let mut h = String::from("##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
        let mut filters: Vec<(&str, &str)> = ALL_FILTERS.iter().map(|(_, n, d)| (*n, *d)).collect();
        filters.sort();
        for (name, desc) in filters {
            h.push_str(&format!("##FILTER=<ID={name},Description=\"{desc}\">\n"));
        }
        h.push_str(HEADER_FIELDS);
        for (name, len) in contigs {
            h.push_str(&format!("##contig=<ID={name},length={len}>\n"));
        }
        h.push_str(&format!("##sageVersion={SAGE_VERSION}\n"));
        h.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for s in reference_ids.iter().chain(tumor_ids.iter()) {
            h.push('\t');
            h.push_str(s);
        }
        h.push('\n');
        out.write_all(h.as_bytes())?;
        Ok(VcfWriter { path: path.to_string(), out, pending: Vec::new() })
    }

    fn write_pending(&mut self) -> Result<()> {
        for (_, line) in self.pending.drain(..) {
            self.out.write_all(line.as_bytes())?;
            self.out.write_all(b"\n")?;
        }
        Ok(())
    }

    /// Add one region's records (in order); the previous region is written once its successor is known.
    pub fn add_region(&mut self, mut records: Vec<(i32, String)>) -> Result<()> {
        // checkLaterVariants
        if !self.pending.is_empty() && !records.is_empty() {
            while !records.is_empty() {
                let last = self.pending.len() - 1;
                if records[0].0 >= self.pending[last].0 {
                    break;
                }
                let rec = records.remove(0);
                let mut index = last as i64;
                while index >= 0 {
                    if self.pending[index as usize].0 <= rec.0 {
                        break;
                    }
                    index -= 1;
                }
                self.pending.insert((index + 1) as usize, rec);
            }
        }
        self.write_pending()?;
        self.pending = records;
        Ok(())
    }

    pub fn flush_chromosome(&mut self) -> Result<()> {
        self.write_pending()
    }

    pub fn close(mut self) -> Result<()> {
        self.write_pending()?;
        self.out.flush()?;
        let path = self.path.clone();
        drop(self.out);
        let c_path = std::ffi::CString::new(path.clone())?;
        // tabix index, as htsjdk writes alongside a .vcf.gz
        let ret = unsafe { rust_htslib::htslib::tbx_index_build(c_path.as_ptr(), 0, &rust_htslib::htslib::tbx_conf_vcf) };
        if ret != 0 {
            anyhow::bail!("failed to build tabix index for {path}");
        }
        Ok(())
    }
}
