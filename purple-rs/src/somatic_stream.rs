//! SomaticPurityEnrichment and SomaticStream: the somatic variants annotated with copy number, VAF,
//! biallelic probability, kataegis ids, subclonal likelihood and genotypes, counted for TML, TMB and
//! MSI, checked for reportability, and written to purple.somatic.vcf.gz; the clonality peak model and
//! the copy-number histogram for the plots.

use anyhow::Result;
use std::collections::{BTreeMap, HashMap};

use hmf_common::commons_math::Poisson;
use hmf_common::variant_context::{self as vcf, EncodeOptions, Value, VcfHeader, VcfWriter};

use crate::config::PurpleConfig;
use crate::constants::*;
use crate::copynumber::PurpleCopyNumber;
use crate::drivers::{reportable_transcript_list, DriverCatalog, DriverFinder, DriverSource, VariantTypeCounts};
use crate::fitting::PurityAdjuster;
use crate::fittingsnv::{generate_model_peaks, PeakModelData, SubclonalLikelihood};
use crate::gene::GeneCopyNumber;
use crate::reference::ReferenceData;
use crate::region::ObservedRegion;
use crate::types::*;
use crate::variant::Variant;

/// a region selector over a chromosome-sorted list of (chromosome, start, end) items: the item containing a position
pub struct RegionLookup<T> {
    by_chr: HashMap<String, Vec<(i32, i32, T)>>,
}

impl<T: Clone> RegionLookup<T> {
    pub fn new(items: impl IntoIterator<Item = (String, i32, i32, T)>) -> RegionLookup<T> {
        let mut by_chr: HashMap<String, Vec<(i32, i32, T)>> = HashMap::new();
        for (c, s, e, t) in items {
            if Chromosome::contains(&c) {
                by_chr.entry(chr(&c).short_name().to_string()).or_default().push((s, e, t));
            }
        }
        for list in by_chr.values_mut() {
            list.sort_by_key(|x| x.0);
        }
        RegionLookup { by_chr }
    }
    pub fn get(&self, chromosome: &str, position: i32) -> Option<&T> {
        let list = self.by_chr.get(chr(chromosome).short_name())?;
        let i = list.partition_point(|x| x.0 <= position);
        if i == 0 {
            return None;
        }
        let (_, e, t) = &list[i - 1];
        if position <= *e { Some(t) } else { None }
    }
}

// ---------------------------------------------------------------------------------------------
// SomaticPurityEnrichment

/// SomaticPurityEnrichment.calculateBiallelic
pub fn calculate_biallelic(cn: &PurpleCopyNumber, allele_read_count: i32, variant_copy_number: f64) -> f64 {
    if allele_read_count == 0 {
        return 0.0;
    }
    let minor = cn.minor_allele_copy_number();
    let copy_number = cn.average_tumor_copy_number;
    let p_loh = 1.0 - 1.0 / (1.0 + jexp(-(BIALLELIC_LOH_GROWTH_RATE as f64) * (minor - 0.5)));
    let p_no_loh = 1.0 - p_loh;
    let threshold = jmin(copy_number - 0.5, jmax(BIALLELIC_THRESHOLD_PARAMETER_I, copy_number - BIALLELIC_THRESHOLD_PARAMETER_II));
    let read_count_at_threshold = (threshold / variant_copy_number) * allele_read_count as f64;
    let p_assume_loh = 1.0 - Poisson::new(allele_read_count as f64).cumulative_probability(read_count_at_threshold.floor() as i32);
    let mut p_assume_no_loh = jmax(p_loh, BIALLELIC_LOH_BASE_ERROR_RATE) / ((1.0 - p_assume_loh) + jmax(p_loh, BIALLELIC_LOH_BASE_ERROR_RATE));
    if p_assume_no_loh.is_nan() {
        p_assume_no_loh = 0.0;
    }
    p_loh * p_assume_loh + p_no_loh * p_assume_no_loh
}

/// SomaticPurityEnrichment.processVariant for every variant
pub fn purity_enrich(variants: &mut [Variant], adjuster: &PurityAdjuster, copy_numbers: &[PurpleCopyNumber], fitted: &[ObservedRegion]) {
    let cn_lookup = RegionLookup::new(copy_numbers.iter().enumerate().map(|(i, c)| (c.chromosome.clone(), c.start, c.end, i)));
    let region_lookup = RegionLookup::new(fitted.iter().map(|r| (r.chromosome.clone(), r.start, r.end, r.germline_status)));
    for v in variants.iter_mut() {
        if !Chromosome::contains(v.chromosome()) {
            continue;
        }
        let status = region_lookup.get(v.chromosome(), v.position()).copied().unwrap_or(GermlineStatus::Unknown);
        v.rec.info.put("PURPLE_GERMLINE", Value::Str(status.name().to_string()));
        if v.tumor_depth.is_none() {
            continue;
        }
        let Some(&ci) = cn_lookup.get(v.chromosome(), v.position()) else { continue };
        let cn = &copy_numbers[ci];
        let copy_number = cn.average_tumor_copy_number;
        let vaf = adjuster.purity_adjusted_vaf_het(&cn.chromosome, jmax(0.001, copy_number), v.allele_frequency(), status == GermlineStatus::HetDeletion);
        let vcn = jmax(0.0, vaf * copy_number);
        v.rec.info.put("PURPLE_VCN", Value::Double(vcn));
        v.rec.info.put("PURPLE_CN", Value::Double(copy_number));
        v.rec.info.put("PURPLE_AF", Value::Str(hmf_common::fmt::java_fixed(vaf, 4)));
        v.rec.info.put("PURPLE_MACN", Value::Double(cn.minor_allele_copy_number()));
        let p = calculate_biallelic(cn, v.allele_read_count(), vcn);
        v.rec.info.put("BIALLELIC_PROB", Value::Str(hmf_common::fmt::java_fixed(p, 4)));
        v.rec.info.put("BIALLELIC", Value::Bool(p >= 0.5));
    }
}

// ---------------------------------------------------------------------------------------------
// kataegis

fn kataegis_forward(v: &Variant) -> bool {
    let alt = v.rec.alts.iter().any(|a| a == "T" || a == "G");
    let tri = v.rec.info.as_string("TNC", "");
    v.is_pass() && tri.starts_with("TC") && alt
}

fn kataegis_reverse(v: &Variant) -> bool {
    let alt = v.rec.alts.iter().any(|a| a == "C" || a == "A");
    let tri = v.rec.info.as_string("TNC", "");
    v.is_pass() && tri.ends_with("GA") && alt
}

/// the id source of the kataegis queues: Java's shared AtomicInteger, or placeholders numbered later
pub struct KataegisIds {
    pub next: i32,
}

impl KataegisIds {
    fn take(&mut self) -> i32 {
        self.next += 1;
        self.next
    }
}

/// KataegisQueue over variant indices
struct KataegisQueue {
    prefix: &'static str,
    candidate: fn(&Variant) -> bool,
    buffer: std::collections::VecDeque<usize>,
    current_id: i32,
}

const MAX_ABS_DISTANCE: i64 = 2000;
const MAX_AVG_DISTANCE: i64 = 1000;
const MIN_COUNT: i32 = 3;

/// (id assigned to a cluster member) recorded as (variant index, prefix, id)
type Assignment = (usize, &'static str, i32);

impl KataegisQueue {
    fn new(prefix: &'static str, forward: bool, ids: &mut KataegisIds) -> KataegisQueue {
        Self::with_predicate(prefix, if forward { kataegis_forward } else { kataegis_reverse }, ids)
    }

    fn with_predicate(prefix: &'static str, candidate: fn(&Variant) -> bool, ids: &mut KataegisIds) -> KataegisQueue {
        let current_id = ids.take();
        KataegisQueue { prefix, candidate, buffer: std::collections::VecDeque::new(), current_id }
    }

    fn candidate(&self, v: &Variant) -> bool {
        (self.candidate)(v)
    }

    /// KataegisQueue.processVariant; variants leaving the queue go to `consumer` at once, as Java's
    /// Consumer does
    fn process(&mut self, index: usize, variants: &[Variant], ids: &mut KataegisIds, out: &mut Vec<Assignment>, consumer: &mut Option<&mut KataegisQueue>) {
        if let Some(&prev) = self.buffer.back() {
            let (p, c) = (&variants[prev], &variants[index]);
            if c.position() as i64 - p.position() as i64 > MAX_ABS_DISTANCE || p.chromosome() != c.chromosome() {
                self.flush(variants, ids, out, consumer);
            }
        }
        self.buffer.push_back(index);
    }

    fn flush(&mut self, variants: &[Variant], ids: &mut KataegisIds, out: &mut Vec<Assignment>, consumer: &mut Option<&mut KataegisQueue>) {
        while !self.buffer.is_empty() {
            self.process_first(variants, ids, out, consumer);
        }
    }

    fn pass_on(index: usize, variants: &[Variant], ids: &mut KataegisIds, out: &mut Vec<Assignment>, consumer: &mut Option<&mut KataegisQueue>) {
        if let Some(c) = consumer.as_deref_mut() {
            c.process(index, variants, ids, out, &mut None);
        }
    }

    fn process_first(&mut self, variants: &[Variant], ids: &mut KataegisIds, out: &mut Vec<Assignment>, consumer: &mut Option<&mut KataegisQueue>) {
        let Some(&first) = self.buffer.front() else { return };
        if !self.candidate(&variants[first]) {
            let v = self.buffer.pop_front().unwrap();
            Self::pass_on(v, variants, ids, out, consumer);
            return;
        }
        // longestViableWindow
        let first_pos = variants[first].position() as i64;
        let (mut w_end, mut w_count) = (first_pos, 0i64);
        let mut best: Option<(i64, i64)> = None;
        for &i in &self.buffer {
            let v = &variants[i];
            if v.position() as i64 - w_end > MAX_ABS_DISTANCE {
                break;
            }
            if self.candidate(v) {
                w_count += 1;
                w_end = v.position() as i64;
            }
            let avg = if w_count <= 1 { 0 } else { hmf_common::java_math::java_round((w_end - first_pos) as f64 / (w_count - 1) as f64) };
            if w_count >= MIN_COUNT as i64 && avg <= MAX_AVG_DISTANCE {
                best = Some((w_end, w_count));
            }
        }
        // the result window keeps the end and count; its start is set to its end (KataegisWindow's copy constructor)
        let (window_end, viable) = match best {
            Some((end, count)) => (end, count >= MIN_COUNT as i64),
            None => (first_pos, false),
        };
        if viable {
            self.current_id = ids.take();
        }
        while let Some(&peek) = self.buffer.front() {
            if variants[peek].position() as i64 > window_end {
                return;
            }
            if viable && self.candidate(&variants[peek]) {
                out.push((peek, self.prefix, self.current_id));
            }
            let v = self.buffer.pop_front().unwrap();
            Self::pass_on(v, variants, ids, out, consumer);
        }
    }
}

/// KataegisEnrichment over one enricher's variants in order: (variant, prefix, id) assignments
fn kataegis_enricher(order: &[usize], variants: &[Variant], fwd: &mut KataegisQueue, rev: &mut KataegisQueue, ids: &mut KataegisIds, out: &mut Vec<Assignment>) {
    for &i in order {
        fwd.process(i, variants, ids, out, &mut Some(rev));
    }
    fwd.flush(variants, ids, out, &mut Some(rev));
    rev.flush(variants, ids, out, &mut None);
}

/// The KT ids. Java (PU1) runs `threads` enrichers in parallel, each with its own queues (each queue
/// takes an id when it is built) and the chromosomes dealt round-robin, all drawing from one counter, so
/// the ids depend on thread timing. Java mode reproduces the enrichers run one after the other in task
/// order (Java's result at -threads 1, and one of its possible results at more threads); the default
/// numbers the clusters 1, 2, ... in genome order of their first variant.
pub fn assign_kataegis(variants: &mut [Variant], enriched: &[usize], threads: i32) {
    let java = crate::compat::java("PU1");
    let mut out: Vec<Assignment> = Vec::new();
    let mut ids = KataegisIds { next: 0 };
    if java && threads > 1 {
        let n = threads as usize;
        let mut tasks: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut task = 0;
        let mut current = variants.first().map(|v| v.chromosome().to_string()).unwrap_or_default();
        for (i, v) in variants.iter().enumerate() {
            if current != v.chromosome() {
                current = v.chromosome().to_string();
                task += 1;
                if task >= n {
                    task = 0;
                }
            }
            tasks[task].push(i);
        }
        let mut queues: Vec<(KataegisQueue, KataegisQueue)> = Vec::new();
        for _ in 0..n {
            let rev = KataegisQueue::new("REV", false, &mut ids);
            let fwd = KataegisQueue::new("FWD", true, &mut ids);
            queues.push((fwd, rev));
        }
        let enriched_set: std::collections::HashSet<usize> = enriched.iter().copied().collect();
        for (t, (fwd, rev)) in queues.iter_mut().enumerate() {
            let order: Vec<usize> = tasks[t].iter().copied().filter(|i| enriched_set.contains(i)).collect();
            kataegis_enricher(&order, variants, fwd, rev, &mut ids, &mut out);
        }
    } else {
        let mut rev = KataegisQueue::new("REV", false, &mut ids);
        let mut fwd = KataegisQueue::new("FWD", true, &mut ids);
        kataegis_enricher(enriched, variants, &mut fwd, &mut rev, &mut ids, &mut out);
    }
    if !java {
        // renumber: clusters by (genome order of their first member)
        let mut first: BTreeMap<(&'static str, i32), usize> = BTreeMap::new();
        for (i, p, id) in &out {
            let e = first.entry((*p, *id)).or_insert(*i);
            if *i < *e {
                *e = *i;
            }
        }
        let mut clusters: Vec<((&'static str, i32), usize)> = first.into_iter().collect();
        clusters.sort_by_key(|(_, i)| *i);
        let renumber: HashMap<(&'static str, i32), i32> = clusters.iter().enumerate().map(|(n, (k, _))| (*k, n as i32 + 1)).collect();
        for a in out.iter_mut() {
            a.2 = renumber[&(a.1, a.2)];
        }
    }
    for (i, p, id) in out {
        variants[i].rec.info.put("KT", Value::Str(format!("{p}_{id}")));
    }
}

// ---------------------------------------------------------------------------------------------
// TML, MSI

/// MicrosatelliteIndels.repeatContextIsRelevant
fn repeat_context_relevant(count: i32, seq_len: usize) -> bool {
    let long = (2..=4).contains(&seq_len) && count >= 4;
    let short = seq_len == 1 && count >= 5;
    long | short
}

pub struct SomaticStreamResult {
    pub tml: f64,
    pub tmb: f64,
    pub msi_indels_per_mb: f64,
    pub msi_indel_count: i32,
    pub burden: i32,
    pub enabled: bool,
    pub peak_model: Vec<PeakModelData>,
    pub reported_genes: std::collections::HashSet<String>,
    /// indices of the variants kept for the Circos plots
    pub plotting_variants: Vec<usize>,
    /// the SNVs of the rainfall plot (RChartData's somatic_data.tsv)
    pub rainfall: Vec<crate::plots::RainfallPoint>,
    pub variant_type_counts: VariantTypeCounts,
    pub biallelic_counts: VariantTypeCounts,
}

impl SomaticStreamResult {
    pub fn empty() -> SomaticStreamResult {
        SomaticStreamResult {
            tml: 0.0,
            tmb: 0.0,
            msi_indels_per_mb: 0.0,
            msi_indel_count: 0,
            burden: 0,
            enabled: false,
            peak_model: Vec::new(),
            reported_genes: Default::default(),
            plotting_variants: Vec::new(),
            rainfall: Vec::new(),
            variant_type_counts: VariantTypeCounts::default(),
            biallelic_counts: VariantTypeCounts::default(),
        }
    }
    pub fn tumor_mutational_load(&self) -> i32 {
        jround(self.tml) as i32
    }
    pub fn microsatellite_status(&self) -> MicrosatelliteStatus {
        if self.enabled { MicrosatelliteStatus::from_indels_per_mb(self.msi_indels_per_mb) } else { MicrosatelliteStatus::Unknown }
    }
    pub fn tmb_status(&self) -> TumorMutationalStatus {
        if self.enabled { TumorMutationalStatus::from_burden_per_mb(self.tmb) } else { TumorMutationalStatus::Unknown }
    }
    pub fn tml_status(&self) -> TumorMutationalStatus {
        if self.enabled { TumorMutationalStatus::from_load(self.tml) } else { TumorMutationalStatus::Unknown }
    }
}

fn has_phased_effect(v: &Variant) -> bool {
    let e = &v.impact.canonical_effect;
    e.contains("phased_inframe_insertion") || e.contains("phased_inframe_deletion") || e.contains("phased_missense")
}

pub struct SomaticStream<'a> {
    pub config: &'a PurpleConfig,
    pub reference: &'a ReferenceData,
    pub onco: DriverFinder<'a>,
    pub tsg: DriverFinder<'a>,
}

/// the somatic VCF header: the input's plus KT, SUBCL, the hotspot flags and the PURPLE fields
pub fn somatic_header(input: &VcfHeader, version: &str, tumor_only: bool) -> VcfHeader {
    let mut h = input.clone();
    h.add_line(vcf::info_line("KT", "1", "String", "Forward/reverse kataegis id"));
    h.add_line(vcf::info_line("SUBCL", "1", "Float", "Non-zero subclonal likelihood"));
    h.add_line(vcf::info_line("HOTSPOT", "0", "Flag", "Site is at a known hotspot location"));
    h.add_line(vcf::info_line("NEAR_HOTSPOT", "0", "Flag", "Variant within 5 bases of hotspot"));
    h.add_line(vcf::meta_line("purpleVersion", version));
    h.add_line(vcf::info_line("PURPLE_AF", "1", "Float", "Purity adjusted variant allelic frequency"));
    h.add_line(vcf::info_line("PURPLE_CN", "1", "Float", "Purity adjusted copy number surrounding variant location"));
    h.add_line(vcf::info_line("PURPLE_VCN", "1", "Float", "Purity adjusted variant copy number"));
    h.add_line(vcf::info_line("PURPLE_MACN", "1", "Float", "Purity adjusted minor allele ploidy surrounding variant location"));
    h.add_line(vcf::info_line("PURPLE_GERMLINE", "1", "String", "Germline classification surrounding variant location"));
    h.add_line(vcf::info_line("BIALLELIC", "0", "Flag", "If the probability this variant is biallelic equals or is greater than 0.50"));
    h.add_line(vcf::info_line("BIALLELIC_PROB", "1", "Float", "Biallelic probability"));
    h.add_line(vcf::info_line("REPORTED", "0", "Flag", "Variant is reported in the driver catalog"));
    h.add_line(vcf::info_line("REPORTABLE_TRANSCRIPTS", "1", "String", "List of reportable transcript when non-canonical are reportable"));
    if tumor_only {
        h.add_line(vcf::info_line("VAF_DIS_MIN", "2", "Float", "Panel germline VAF distance, somatic VAF distance"));
        h.add_line(vcf::info_line("SOM_LH", "1", "String", "Panel somatic likelihood [HIGH,MEDIUM,LOW]"));
    }
    h
}

impl<'a> SomaticStream<'a> {
    pub fn new(config: &'a PurpleConfig, reference: &'a ReferenceData) -> SomaticStream<'a> {
        SomaticStream {
            config,
            reference,
            onco: DriverFinder::new(&reference.driver_genes, DriverCategory::Onco),
            tsg: DriverFinder::new(&reference.driver_genes, DriverCategory::Tsg),
        }
    }

    /// SomaticStream.processAndWrite
    #[allow(clippy::too_many_arguments)]
    pub fn process_and_write(
        &mut self,
        cache: &mut crate::somatic::SomaticVariantCache,
        adjuster: Option<&PurityAdjuster>,
        msi_prediction: Option<f64>,
        output_vcf: &str,
        hist_path: &str,
        collect_rainfall: bool,
    ) -> Result<SomaticStreamResult> {
        let mut result = SomaticStreamResult::empty();
        let Some(input_header) = cache.header.clone() else { return Ok(result) };
        result.enabled = true;
        crate::debug!("modelling somatic peaks");
        result.peak_model = generate_model_peaks(&cache.variants);
        let config = self.config;
        let tumor_only = config.tumor_only_mode();
        let header = somatic_header(&input_header, &config.version, tumor_only);

        // enrichment: kataegis, subclonal likelihood, genotypes
        let variants = &mut cache.variants;
        let enriched: Vec<usize> = (0..variants.len()).filter(|&i| !(tumor_only && variants[i].is_filtered() && !config.write_all_somatics)).collect();
        assign_kataegis(variants, &enriched, config.threads);
        let subclonal = SubclonalLikelihood::new(CLONALITY_BIN_WIDTH, &result.peak_model);
        for &i in &enriched {
            let v = &mut variants[i];
            let l = jround(subclonal.subclonal_likelihood(v.copy_number()) * 1000.0) as f64 / 1000.0;
            v.rec.info.put("SUBCL", Value::Double(l));
            genotype_enrich(v, &input_header, config.reference_id.as_deref(), &config.tumor_id);
        }

        // TML, MSI, reportability, charts
        let targets = &self.reference.target_regions;
        let (mut load, mut burden, mut unclear, mut msi_count) = (0i32, 0i32, 0i32, 0i32);
        let mut histogram: HashMap<String, i32> = HashMap::new();
        // RChartData: the SNVs of the rainfall plot, at most 100,000 (PL4: Java counts every PASS variant)
        let mut rainfall_seen = 0usize;
        let (snp_mod, indel_mod) = {
            let m = |c: i32| if c <= 25000 { 1 } else { c / 25000 };
            (m(cache.snp_count), m(cache.indel_count))
        };
        let (mut snp_count, mut indel_count) = (0i32, 0i32);
        for i in 0..variants.len() {
            if !Chromosome::contains(variants[i].chromosome()) {
                continue;
            }
            if let Some(adj) = adjuster {
                if tumor_only && (variants[i].is_pass() || config.write_all_somatics) {
                    germline_likelihood(&mut variants[i], &input_header, cache.tumor_ordinal, adj.purity);
                }
            }
            if !variants[i].is_pass() {
                continue;
            }
            let v = &variants[i];
            // TumorMutationalLoad
            if v.allele_frequency() >= TUMOR_MUT_LOAD_MIN_VAF {
                if targets.has_target_regions() {
                    if targets.in_target_regions(v.chromosome(), v.position())
                        && !v.is_hotspot()
                        && v.vtype == VariantType::Snp
                        && v.impact.worst_coding_effect == CodingEffect::Missense
                        && !TARGETED_TMB_GENE_EXCLUSIONS.contains(&v.impact.gene_name.as_str())
                        && v.rec.info.as_double("GND_FREQ", 0.0) <= 0.0
                    {
                        let mut is_unclear = false;
                        let mut skip = false;
                        if tumor_only {
                            let lh = if v.rec.info.has("SOM_LH") { v.rec.info.as_string("SOM_LH", "") } else { "LOW".into() };
                            if lh != "HIGH" {
                                if lh == "MEDIUM" { is_unclear = true } else { skip = true }
                            }
                        }
                        if !skip {
                            if is_unclear { unclear += 1 } else { load += 1 }
                        }
                    }
                } else {
                    burden += 1;
                    if v.impact.worst_coding_effect == CodingEffect::Missense {
                        load += 1;
                    }
                }
            }
            // MicrosatelliteIndels
            if v.vtype == VariantType::Indel && v.ref_bases().len() < 50 && v.alt().len() < 50 && v.allele_frequency() >= TUMOR_MSI_LOAD_MIN_VAF {
                let rc = v.rec.info.as_int("REP_C", 0);
                let seq = v.rec.info.as_string("REP_S", "");
                if repeat_context_relevant(rc, seq.len()) {
                    msi_count += 1;
                }
            }
            // checkReportability
            let reportable = {
                let v = &variants[i];
                let mut r = false;
                if !v.is_filtered() {
                    result.variant_type_counts.increment(v.vtype);
                    if v.biallelic() {
                        result.biallelic_counts.increment(v.vtype);
                    }
                    r = self.onco.check_variant(i, v) || self.tsg.check_variant(i, v);
                }
                r
            };
            if reportable {
                let v = &mut variants[i];
                v.rec.info.put("REPORTED", Value::Bool(true));
                result.reported_genes.insert(v.gene().to_string());
                if let Some(t) = reportable_transcript_list(v.vtype, &v.impact) {
                    v.rec.info.put("REPORTABLE_TRANSCRIPTS", Value::Str(t));
                }
            }
            // RChartData.processVariant
            let v = &variants[i];
            let cn = v.rec.info.as_double("PURPLE_CN", 0.0);
            let vcn = v.rec.info.as_double("PURPLE_VCN", 0.0);
            let key = format!("{}>{}", jround(vcn / 0.05) as i32, jround(cn / 1.0) as i32);
            *histogram.entry(key).or_insert(0) += 1;
            if collect_rainfall {
                let snv = v.vtype == VariantType::Snp && v.is_pass();
                let counted = if crate::compat::java("PL4") { rainfall_seen < 100_000 } else { result.rainfall.len() < 100_000 };
                if counted {
                    rainfall_seen += 1;
                    if snv {
                        result.rainfall.push(crate::plots::RainfallPoint {
                            chromosome: v.chromosome().to_string(),
                            position: v.position() as i64,
                            mutation: format!("{}>{}", v.ref_bases(), v.alt()),
                            kataegis: v.rec.info.as_string("KT", ""),
                        });
                    }
                }
            }
            // checkChartDownsampling
            if !config.charting.disabled {
                let reported = v.reported();
                if !(!reported && targets.has_target_regions() && !targets.in_target_regions(v.chromosome(), v.position())) {
                    if v.vtype == VariantType::Indel {
                        indel_count += 1;
                        if reported || indel_count % indel_mod == 0 {
                            result.plotting_variants.push(i);
                        }
                    } else {
                        snp_count += 1;
                        if reported || snp_count % snp_mod == 0 {
                            result.plotting_variants.push(i);
                        }
                    }
                }
            }
        }
        self.check_phased_reportable(variants);

        // write the VCF
        let mut writer = VcfWriter::create(output_vcf, &header, EncodeOptions { allow_missing_fields: true })?;
        for v in variants.iter() {
            if !tumor_only || v.is_pass() || config.write_all_somatics {
                writer.write(&v.rec)?;
            }
        }
        writer.close()?;

        // RChartData.write
        let mut lines: Vec<String> = histogram
            .iter()
            .map(|(k, c)| {
                let (a, b) = k.split_once('>').unwrap();
                let a: i32 = a.parse().unwrap();
                let b: i32 = b.parse().unwrap();
                format!("{}\t{}\t{}", hmf_common::fmt::java_fixed(a as f64 * 0.05, 2), hmf_common::fmt::java_fixed(b as f64 * 1.0, 0), c)
            })
            .collect();
        lines.sort();
        let mut text = String::from("variantCopyNumberBucket\tcopyNumberBucket\tcount\n");
        for l in lines {
            text.push_str(&l);
            text.push('\n');
        }
        std::fs::write(hist_path, text)?;

        // calculateVariantLoadValues
        result.tml = if !targets.has_target_regions() {
            load as f64
        } else if targets.coding_bases == 0 {
            0.0
        } else {
            let mut unclear_load = 0.0;
            if unclear > 0 {
                let factor = (targets.coding_bases / targets.coding_base_factor) as f64;
                unclear_load = jpow(unclear as f64, 2.0) / (unclear as f64 + factor);
            }
            (load as f64 + unclear_load) * targets.tml_ratio * CODING_BASES_PER_GENOME / targets.coding_bases as f64
        };
        result.msi_indels_per_mb = if targets.has_target_regions() { msi_prediction.map(|p| jmax(p, 0.0)).unwrap_or(0.0) } else { msi_count as f64 / MB_PER_GENOME };
        result.msi_indel_count = msi_count;
        result.burden = burden;
        result.tmb = if targets.has_target_regions() { result.msi_indels_per_mb + result.tml * targets.tmb_ratio } else { burden as f64 / MB_PER_GENOME };
        crate::info!(
            "load(tml={}) msiIndels({} perMb={}) burden({} perMb={})",
            hmf_common::fmt::java_fixed(result.tml, 4),
            msi_count,
            hmf_common::fmt::java_fixed(result.msi_indels_per_mb, 4),
            hmf_common::fmt::java_fixed(burden as f64, 1),
            hmf_common::fmt::java_fixed(result.tmb, 4)
        );
        Ok(result)
    }

    /// SomaticStream.checkPhasedReportableVariants
    fn check_phased_reportable(&mut self, variants: &mut [Variant]) {
        for i in 0..variants.len() {
            if !variants[i].rec.info.has("REPORTED") || !has_phased_effect(&variants[i]) {
                continue;
            }
            let lps = variants[i].rec.info.as_int_list("LPS", 0);
            if lps.is_empty() {
                continue;
            }
            for back in [true, false] {
                let mut j = i as i64;
                loop {
                    j += if back { -1 } else { 1 };
                    if j < 0 || j as usize >= variants.len() {
                        break;
                    }
                    let ju = j as usize;
                    let next = &variants[ju];
                    // PU8: Java tests the reported variant's phased effect again here instead of the neighbour's
                    let phased = if crate::compat::java("PU8") { has_phased_effect(&variants[i]) } else { has_phased_effect(next) };
                    if !next.is_pass() || next.rec.info.has("REPORTED") || !phased {
                        continue;
                    }
                    if !next.impact.canonical_coding_effect.has_protein_impact() {
                        continue;
                    }
                    let next_lps = next.rec.info.as_int_list("LPS", 0);
                    if next_lps.is_empty() || !next_lps.iter().any(|x| lps.contains(x)) {
                        break;
                    }
                    variants[ju].rec.info.put("REPORTED", Value::Bool(true));
                    crate::debug!("var({}) setting reported due to inframe-phasing with other({})", variants[ju].describe(), variants[i].describe());
                    if self.onco.driver_variants.contains(&i) {
                        self.onco.driver_variants.push(ju);
                    }
                    if self.tsg.driver_variants.contains(&i) {
                        self.tsg.driver_variants.push(ju);
                    }
                }
            }
        }
    }

    /// SomaticStream.buildDrivers (with the target-region count override)
    pub fn build_drivers(
        &self,
        variants: &[Variant],
        result: &SomaticStreamResult,
        gene_copy_numbers: &HashMap<String, Vec<GeneCopyNumber>>,
        sources: &mut Vec<DriverSource>,
    ) -> Vec<DriverCatalog> {
        let mut counts = result.variant_type_counts.clone();
        let mut biallelic = result.biallelic_counts.clone();
        let targets = &self.reference.target_regions;
        if targets.has_target_regions() {
            let snv = jround(result.tml * targets.tmb_ratio * MB_PER_GENOME) as i32;
            let indel = jround(result.msi_indels_per_mb * MB_PER_GENOME) as i32;
            counts.snp = Some(snv);
            counts.indel = Some(indel);
            biallelic.snp = Some((snv as f64 * BIALLELIC_ASSUMED_FRACTION) as i32);
            biallelic.indel = Some((indel as f64 * BIALLELIC_ASSUMED_FRACTION) as i32);
            crate::debug!("target-regions inferred driver variant counts snv({snv}) indel({indel})");
        }
        let mut out = self.onco.find_drivers(variants, gene_copy_numbers, &counts, &biallelic, sources);
        out.extend(self.tsg.find_drivers(variants, gene_copy_numbers, &counts, &biallelic, sources));
        out
    }
}

/// SomaticGenotypeEnrichment.processVariant
fn genotype_enrich(v: &mut Variant, header: &VcfHeader, reference: Option<&str>, tumor: &str) {
    if v.rec.alts.is_empty() {
        return;
    }
    let biallelic = v.biallelic();
    let genotypes = v.rec.genotypes(header);
    let mut updated = Vec::new();
    if let Some(r) = reference.filter(|r| !r.is_empty()) {
        if let Some(g) = genotypes.iter().find(|g| g.sample == r) {
            let mut g = g.clone();
            g.alleles = vec![Some(0), Some(0)];
            updated.push(g);
        }
    }
    let t = genotypes.iter().find(|g| g.sample == tumor).cloned().unwrap_or_else(|| panic!("NullPointerException: no genotype for {tumor}"));
    let mut t = t;
    t.alleles = if biallelic { vec![Some(1), Some(1)] } else { vec![Some(0), Some(1)] };
    updated.push(t);
    *genotypes = updated;
    if v.rec.info.has("REPORTED") {
        v.rec.info.remove("REPORTED");
    }
    v.refresh();
}

/// SomaticGermlineLikelihood.processVariant (tumor-only mode)
fn germline_likelihood(v: &mut Variant, header: &VcfHeader, tumor_ordinal: i32, purity: f64) {
    let raw_af = {
        let gs = v.rec.genotypes(header);
        gs.get(tumor_ordinal as usize).map(|g| g.attr_double("AF", 0.0)).unwrap_or(0.0)
    };
    let segment_cn = v.rec.info.as_double("PURPLE_CN", 0.0);
    let tumor_minor = v.rec.info.as_double("PURPLE_MACN", 0.0);
    let tumor_major = segment_cn - tumor_minor;
    let ref_purity = 1.0 - purity;
    let denom = 2.0 * ref_purity + segment_cn * purity;
    let g_major = (ref_purity + tumor_major * purity) / denom;
    let g_minor = (ref_purity + tumor_minor * purity) / denom;
    let (gmaj, gmin, ghom) = (g_major - raw_af, g_minor - raw_af, 1.0 - raw_af);
    let mut gdiff = if gmaj.abs() < gmin.abs() { gmaj } else { gmin };
    gdiff = if ghom.abs() < gdiff.abs() { ghom } else { gdiff };
    let s1 = purity / denom;
    let smaj = (tumor_major * purity) / denom - raw_af;
    let smin = (tumor_minor * purity) / denom - raw_af;
    let shom = s1 - raw_af;
    let mut sdiff = if smaj.abs() < smin.abs() { smaj } else { smin };
    sdiff = if shom.abs() < sdiff.abs() { shom } else { sdiff };
    if raw_af < s1 {
        sdiff = 0.0;
    }
    v.rec.info.put("VAF_DIS_MIN", Value::DoubleList(vec![gdiff, sdiff]));
    let diff = sdiff.abs() - gdiff.abs();
    let lh = if v.is_hotspot() || diff < PANEL_SOMATIC_LIKELIHOOD_DIFF_HIGH {
        "HIGH"
    } else if diff > PANEL_SOMATIC_LIKELIHOOD_DIFF_LOW {
        "LOW"
    } else {
        "MEDIUM"
    };
    v.rec.info.put("SOM_LH", Value::Str(lh.into()));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::variant;

    // KataegisQueueTest
    fn kataegis(input: &[Variant]) -> Vec<Option<String>> {
        let mut ids = KataegisIds { next: 0 };
        let mut q = KataegisQueue::with_predicate("TST", |v| v.rec.alts[0] == "T", &mut ids);
        let mut out = Vec::new();
        for i in 0..input.len() {
            q.process(i, input, &mut ids, &mut out, &mut None);
        }
        q.flush(input, &mut ids, &mut out, &mut None);
        let mut kt = vec![None; input.len()];
        for (i, prefix, id) in out {
            kt[i] = Some(format!("{prefix}_{id}"));
        }
        kt
    }

    fn create(chromosome: &str, pos: i32, kataegis: bool) -> Variant {
        variant(chromosome, pos, "C", if kataegis { "T" } else { "A" })
    }

    #[test]
    fn expected_behaviour() {
        let input: Vec<Variant> = (0..MIN_COUNT).map(|i| create("1", 100 + i, true)).collect();
        let kt = kataegis(&input);
        assert_eq!(kt.len(), MIN_COUNT as usize);
        assert_eq!(kt[0].as_deref(), Some("TST_2"));
    }

    #[test]
    fn change_contig() {
        let input: Vec<Variant> = (0..MIN_COUNT).map(|i| create(if i % 2 == 0 { "1" } else { "2" }, 100 + i, true)).collect();
        assert_eq!(kataegis(&input)[0], None);
    }

    #[test]
    fn average_distance_clustered() {
        let input: Vec<Variant> = [100, 102, 103, 104, 105, 2100].iter().map(|&p| create("1", p, true)).collect();
        assert_eq!(kataegis(&input)[0].as_deref(), Some("TST_2"));
    }

    #[test]
    fn average_distance_spread() {
        let input: Vec<Variant> = [1101, 2102, 3103, 4104, 5105, 6103].iter().map(|&p| create("1", p, true)).collect();
        assert_eq!(kataegis(&input)[0].as_deref(), Some("TST_2"));
    }

    #[test]
    fn max_depth_acceptable() {
        let mut input: Vec<Variant> = (0..MIN_COUNT).map(|i| create("1", 100 + i, true)).collect();
        let i = MIN_COUNT - 1;
        input.push(create("1", 100 + i + MAX_ABS_DISTANCE as i32, true));
        input.push(create("1", 100 + i + 2 * MAX_ABS_DISTANCE as i32 + 1, true));
        let kt = kataegis(&input);
        assert_eq!(kt.len(), MIN_COUNT as usize + 2);
        assert_eq!(kt[i as usize + 1].as_deref(), Some("TST_2"));
        assert_eq!(kt[i as usize + 2], None);
    }
}
