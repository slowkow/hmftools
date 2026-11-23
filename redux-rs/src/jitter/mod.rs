//! Microsatellite jitter analysis and the MSI prediction (Java packages redux.jitter and redux.ms_model).
//!
//! Java keeps one synchronized object with boxed hash maps for each of the 2.4 million sites and finds the sites
//! of a read by querying an interval tree, which builds lists, pairs and a stream for every read. Here the sites
//! of a chromosome are arrays sorted by start, each worker follows the sorted reads with a cursor into them, and
//! counts go into a small dense window owned by the worker that is added to the shared totals a few thousand
//! sites at a time. The result does not depend on the thread count: every count is a sum over reads.

mod fit;
mod ms_model;
mod read;
mod table;

use crate::config::{Config, SequencingType};
use anyhow::{anyhow, bail, Context, Result};
use fit::SearchGrid;
use read::{ConsensusType, SiteRef};
use rust_htslib::bam;
use std::collections::{BTreeMap, HashMap};
use std::io::{Read as _, Write};
use std::sync::Mutex;
use table::{JitterCountsTable, JitterModelParams};

// JitterConstants
const MIN_PASSING_SITE_READS: i32 = 20;
const ALT_COUNT_FRACTION_INIT: f64 = 0.3;
const ALT_COUNT_FRACTION_STEP: f64 = -0.05;
const MAX_REJECTED_READ_FRACTION: f64 = 0.2;
/// Doubles.EPSILON
const EPSILON: f64 = 1e-10;

// MsJitterConfig defaults
const DEFAULT_MIN_MAPPING_QUALITY: i32 = 50;
const DEFAULT_NUM_SITES_PER_TYPE: i32 = 5000;
const DEFAULT_MAX_SINGLE_SITE_ALT_CONTRIBUTION: f64 = 0.8;

/// Jitters -DENSE_JITTER..=DENSE_JITTER are counted in arrays; anything further out goes to a map.
const DENSE_JITTER: i32 = 5;
const DENSE_SLOTS: usize = (2 * DENSE_JITTER + 1) as usize;
/// Values kept for a site and consensus type: read count, rejected count, then the dense jitter counts.
const TYPE_STRIDE: usize = 2 + DENSE_SLOTS;
/// Values kept for a site in a worker's window: passing count, then TYPE_STRIDE for each consensus type.
const SITE_STRIDE: usize = 1 + 3 * TYPE_STRIDE;
/// Sites a worker counts before it adds them to the shared totals.
const MAX_WINDOW_SITES: usize = 4096;

/// The repeat unit groups of the tables: the six named groups of MsJitterAnalyser, then unit lengths 3, 4 and 5.
const UNIT_GROUPS: [&[&str]; 6] =
    [&["A", "T"], &["C", "G"], &["AT", "TA"], &["AC", "CA", "GT", "TG"], &["AG", "GA", "CT", "TC"], &["CG", "GC"]];
const GROUP_COUNT: usize = 9;
const NO_GROUP: u8 = u8::MAX;

/// MicrosatelliteSelector.unitName for the nine chart tables.
fn group_name(group: usize) -> String {
    if group < 6 { UNIT_GROUPS[group].join("/") } else { format!("{}bp repeat", group - 3) }
}

/// The group a repeat unit is counted in; none for a unit of one or two bases that is in no list, or over five.
fn unit_group(unit: &[u8]) -> u8 {
    for (g, units) in UNIT_GROUPS.iter().enumerate() {
        if units.iter().any(|u| u.as_bytes() == unit) {
            return g as u8;
        }
    }
    match unit.len() {
        3..=5 => (unit.len() + 3) as u8,
        _ => NO_GROUP,
    }
}

/// MsJitterConfig
struct JitterConfig {
    sample_id: String,
    output_dir: String,
    sequencing_type: SequencingType,
    uses_duplex_umis: bool,
    apply_bq_filter: bool,
    /// the lowest mapping quality of a read that is analysed
    min_map_quality: i32,
    write_site_file: bool,
    coefficients_file: Option<String>,
    error_rates_file: Option<String>,
}

impl JitterConfig {
    /// MsJitterConfig.consensusTypes
    fn consensus_types(&self) -> Vec<ConsensusType> {
        let mut types = vec![ConsensusType::None];
        match self.sequencing_type {
            SequencingType::Illumina if self.uses_duplex_umis => types.extend([ConsensusType::Single, ConsensusType::Dual]),
            SequencingType::Sbx | SequencingType::Ultima => types.push(ConsensusType::Dual),
            _ => {}
        }
        types
    }
    /// MsJitterConfig.runMsiPrediction
    fn msi_prediction_files(&self) -> Option<(&str, &str)> {
        match (&self.coefficients_file, &self.error_rates_file) {
            (Some(c), Some(e)) if std::path::Path::new(c).exists() && std::path::Path::new(e).exists() => Some((c, e)),
            _ => None,
        }
    }
}

/// The reference microsatellites (MicrosatelliteSite), grouped by chromosome and sorted by start within each.
struct Sites {
    start: Vec<i32>,
    end: Vec<i32>,
    repeat_count: Vec<i32>,
    group: Vec<u8>,
    unit_offset: Vec<u32>,
    unit_len: Vec<u8>,
    units: Vec<u8>,
    /// position in the sites file
    file_index: Vec<u32>,
    /// chromosome name and its range of site indices
    chromosomes: Vec<(String, usize, usize)>,
}

impl Sites {
    fn len(&self) -> usize {
        self.start.len()
    }
    #[inline]
    fn unit(&self, i: usize) -> &[u8] {
        let o = self.unit_offset[i] as usize;
        &self.units[o..o + self.unit_len[i] as usize]
    }
    #[inline]
    fn site_ref(&self, i: usize) -> SiteRef<'_> {
        SiteRef { start: self.start[i], end: self.end[i], unit: self.unit(i), repeat_count: self.repeat_count[i] }
    }
    fn chromosome_range(&self, chromosome: &str) -> Option<(usize, usize)> {
        self.chromosomes.iter().find(|c| c.0 == chromosome).map(|c| (c.1, c.2))
    }
    fn chromosome_of(&self, i: usize) -> &str {
        let c = self.chromosomes.partition_point(|c| c.2 <= i);
        &self.chromosomes[c].0
    }

    /// From (chromosome, start, end, unit) in file order.
    #[cfg(test)]
    fn build(records: Vec<(String, i32, i32, Vec<u8>)>) -> Result<Sites> {
        let mut raw = RawSites::default();
        for (chromosome, start, end, unit) in &records {
            raw.push(chromosome.as_bytes(), *start, *end, unit)?;
        }
        Ok(Sites::from_raw(raw))
    }

    fn from_raw(raw: RawSites) -> Sites {
        let n = raw.start.len();
        let mut order: Vec<(u32, i32, u32)> = (0..n).map(|i| (raw.chromosome[i], raw.start[i], i as u32)).collect();
        order.sort_unstable();
        let mut sites = Sites {
            start: Vec::with_capacity(n),
            end: Vec::with_capacity(n),
            repeat_count: Vec::with_capacity(n),
            group: Vec::with_capacity(n),
            unit_offset: Vec::with_capacity(n),
            unit_len: Vec::with_capacity(n),
            units: raw.units,
            file_index: Vec::with_capacity(n),
            chromosomes: raw.names.into_iter().map(|c| (c, 0, 0)).collect(),
        };
        for (pos, &(c, _, i)) in order.iter().enumerate() {
            let (c, i) = (c as usize, i as usize);
            if sites.chromosomes[c].2 == 0 {
                sites.chromosomes[c].1 = pos;
            }
            sites.chromosomes[c].2 = pos + 1;
            let (start, end, unit_len) = (raw.start[i], raw.end[i], raw.unit_len[i]);
            let unit_offset = raw.unit_offset[i] as usize;
            sites.start.push(start);
            sites.end.push(end);
            // MicrosatelliteSite: RepeatCount = region.baseLength() / unit.length
            sites.repeat_count.push((end - start + 1) / unit_len as i32);
            sites.group.push(unit_group(&sites.units[unit_offset..unit_offset + unit_len as usize]));
            sites.unit_offset.push(unit_offset as u32);
            sites.unit_len.push(unit_len);
            sites.file_index.push(i as u32);
        }
        sites
    }

    /// RefGenomeMicrosatelliteFile.read, then the specific regions filter of MsJitterAnalyser.loadRefGenomeMicrosatellites.
    fn load(filename: &str, config: &Config) -> Result<Sites> {
        let file = std::fs::File::open(filename).with_context(|| format!("failed to load MS indel sites file({filename})"))?;
        let mut data = Vec::new();
        if filename.ends_with(".gz") {
            flate2::read::MultiGzDecoder::new(file).read_to_end(&mut data)?;
        } else {
            std::io::BufReader::new(file).read_to_end(&mut data)?;
        }
        // a line feed at the end of the file does not start another line
        let mut lines = data.strip_suffix(b"\n").unwrap_or(&data).split(|b| *b == b'\n');
        let header = lines.next().filter(|h| !h.is_empty()).ok_or_else(|| anyhow!("empty MS indel sites file({filename})"))?;
        let header = header.strip_suffix(b"\r").unwrap_or(header);
        let col = |name: &str| header.split(|b| *b == b'\t').position(|c| c == name.as_bytes()).ok_or_else(|| anyhow!("missing column {name} in {filename}"));
        let (chr_index, start_index, end_index, unit_index) = (col("chromosome")?, col("start")?, col("end")?, col("unit")?);
        // SpecificRegions.includeRegion: only the regions filter sites, a chromosome list alone does not
        let regions = &config.specific.regions;
        let invalid = || anyhow!("invalid line in MS indel sites file({filename})");
        let parse_int = |v: &[u8]| -> Result<i32> { Ok(std::str::from_utf8(v)?.parse::<i32>()?) };
        let mut raw = RawSites::default();
        for line in lines {
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            let (mut chromosome, mut start, mut end, mut unit): (Option<&[u8]>, Option<&[u8]>, Option<&[u8]>, Option<&[u8]>) = (None, None, None, None);
            for (i, value) in line.split(|b| *b == b'\t').enumerate() {
                if i == chr_index {
                    chromosome = Some(value);
                } else if i == start_index {
                    start = Some(value);
                } else if i == end_index {
                    end = Some(value);
                } else if i == unit_index {
                    unit = Some(value);
                }
            }
            let chromosome = chromosome.ok_or_else(invalid)?;
            let start = parse_int(start.ok_or_else(invalid)?)?;
            let end = parse_int(end.ok_or_else(invalid)?)?;
            let unit = unit.ok_or_else(invalid)?;
            if !regions.is_empty() {
                let name = std::str::from_utf8(chromosome)?;
                if !regions.iter().any(|r| r.overlaps(name, start, end)) {
                    continue;
                }
            }
            raw.push(chromosome, start, end, unit)?;
        }
        Ok(Sites::from_raw(raw))
    }
}

/// The sites in file order, before sorting.
#[derive(Default)]
struct RawSites {
    names: Vec<String>,
    chromosome: Vec<u32>,
    start: Vec<i32>,
    end: Vec<i32>,
    unit_offset: Vec<u32>,
    unit_len: Vec<u8>,
    units: Vec<u8>,
}

impl RawSites {
    fn push(&mut self, chromosome: &[u8], start: i32, end: i32, unit: &[u8]) -> Result<()> {
        // Java divides the region length by the unit length
        if unit.is_empty() || unit.len() > u8::MAX as usize {
            bail!("invalid microsatellite unit at {}:{start}", String::from_utf8_lossy(chromosome));
        }
        // consecutive lines mostly share the chromosome
        let c = match self.chromosome.last() {
            Some(&c) if self.names[c as usize].as_bytes() == chromosome => c,
            _ => match self.names.iter().position(|n| n.as_bytes() == chromosome) {
                Some(c) => c as u32,
                None => {
                    self.names.push(std::str::from_utf8(chromosome)?.to_string());
                    self.names.len() as u32 - 1
                }
            },
        };
        self.chromosome.push(c);
        self.start.push(start);
        self.end.push(end);
        self.unit_offset.push(self.units.len() as u32);
        self.unit_len.push(unit.len() as u8);
        self.units.extend_from_slice(unit);
        Ok(())
    }
}

/// The counts of MicrosatelliteSiteData for every site, summed over the workers. Counts are Java ints.
struct Totals {
    /// mPassingReadCount by site
    passing: Vec<u32>,
    /// by consensus type, TYPE_STRIDE values for each site; empty until a read of the type is counted
    by_type: [Vec<u32>; 3],
    /// (site, consensus type, jitter) -> reads, for jitters outside the dense range
    far_jitter: BTreeMap<(u32, u8, i32), u32>,
    /// (site, consensus type, repeat length in bases) -> reads; kept only for the site file
    repeat_lengths: BTreeMap<(u32, u8, i32), u32>,
}

impl Totals {
    fn new(site_count: usize) -> Self {
        Totals { passing: vec![0; site_count], by_type: [Vec::new(), Vec::new(), Vec::new()], far_jitter: BTreeMap::new(), repeat_lengths: BTreeMap::new() }
    }
    #[inline]
    fn type_values(&self, site: usize, ct: ConsensusType) -> Option<&[u32]> {
        let v = &self.by_type[ct as usize];
        if v.is_empty() { None } else { Some(&v[site * TYPE_STRIDE..(site + 1) * TYPE_STRIDE]) }
    }
    /// MicrosatelliteSiteData.readCountByConsensus
    fn read_count(&self, site: usize, ct: ConsensusType) -> i32 {
        self.type_values(site, ct).map(|v| v[0] as i32).unwrap_or(0)
    }
    /// MicrosatelliteSiteData.readsRejectedByConsensus
    fn rejected_count(&self, site: usize, ct: ConsensusType) -> i32 {
        self.type_values(site, ct).map(|v| v[1] as i32).unwrap_or(0)
    }
    /// MicrosatelliteSiteData.passingJitterCounts: every (jitter, reads) with at least one read.
    fn for_each_jitter(&self, site: usize, ct: ConsensusType, mut f: impl FnMut(i32, i32)) {
        if let Some(v) = self.type_values(site, ct) {
            for (k, &count) in v[2..].iter().enumerate() {
                if count != 0 {
                    f(k as i32 - DENSE_JITTER, count as i32);
                }
            }
        }
        if !self.far_jitter.is_empty() {
            for (&(_, _, jitter), &count) in self.far_jitter.range((site as u32, ct as u8, i32::MIN)..=(site as u32, ct as u8, i32::MAX)) {
                f(jitter, count as i32);
            }
        }
    }
    /// MicrosatelliteSiteData.passingJitterCount
    fn jitter_count(&self, site: usize, ct: ConsensusType, jitter: i32) -> i32 {
        if (-DENSE_JITTER..=DENSE_JITTER).contains(&jitter) {
            self.type_values(site, ct).map(|v| v[2 + (jitter + DENSE_JITTER) as usize] as i32).unwrap_or(0)
        } else if self.far_jitter.is_empty() {
            0
        } else {
            self.far_jitter.get(&(site as u32, ct as u8, jitter)).map(|c| *c as i32).unwrap_or(0)
        }
    }

    /// MicrosatelliteSiteData.isRealVariant with the JitterConstants thresholds.
    #[cfg(test)]
    fn is_real_variant(&self, site: usize) -> bool {
        self.is_real_variant_with(site, crate::compat::java("R11"))
    }

    fn is_real_variant_with(&self, site: usize, java_total: bool) -> bool {
        let passing = self.passing[site] as i32;
        if passing < MIN_PASSING_SITE_READS {
            return true;
        }
        // R11: Java never increments mTotalReadCount, so the rejected fraction is 1 - passing / 0.0 = -Infinity and
        // the test never fires. Corrected: the total is every read of the site, passing or rejected.
        let total_read_count = if java_total {
            0
        } else {
            ConsensusType::ALL.iter().fold(passing, |sum, ct| sum.wrapping_add(self.rejected_count(site, *ct)))
        };
        let fraction_rejected = 1.0 - passing as f64 / total_read_count as f64;
        if fraction_rejected - MAX_REJECTED_READ_FRACTION > -EPSILON {
            return true;
        }
        // reads by jitter over all consensus types, whichever types the run reports
        let mut dense = [0i32; DENSE_SLOTS];
        let mut far: Vec<(i32, i32)> = Vec::new();
        for ct in ConsensusType::ALL {
            self.for_each_jitter(site, ct, |jitter, count| {
                if (-DENSE_JITTER..=DENSE_JITTER).contains(&jitter) {
                    dense[(jitter + DENSE_JITTER) as usize] = dense[(jitter + DENSE_JITTER) as usize].wrapping_add(count);
                } else if let Some(e) = far.iter_mut().find(|e| e.0 == jitter) {
                    e.1 = e.1.wrapping_add(count);
                } else {
                    far.push((jitter, count));
                }
            });
        }
        let over = |repeat_diff: i32, read_count: i32| {
            let fraction_cutoff = fit::java_max(ALT_COUNT_FRACTION_INIT + (repeat_diff.wrapping_abs() - 1) as f64 * ALT_COUNT_FRACTION_STEP, 0.1);
            let count_cutoff = fraction_cutoff * passing as f64;
            read_count as f64 - count_cutoff > EPSILON
        };
        for (k, &count) in dense.iter().enumerate() {
            let jitter = k as i32 - DENSE_JITTER;
            if jitter != 0 && count != 0 && over(jitter, count) {
                return true;
            }
        }
        far.iter().any(|&(jitter, count)| over(jitter, count))
    }
}

/// MsJitterAnalyser: shared by all threads; writes `<sample>.redux.ms_table.tsv.gz`, `<sample>.redux.jitter_params.tsv`
/// and, when the model files are given, the MSI prediction.
pub struct MsJitterAnalyser {
    config: JitterConfig,
    sites: Sites,
    totals: Mutex<Totals>,
}

impl MsJitterAnalyser {
    /// None when `-ref_genome_msi_file` is not given (Java: MsJitterConfig.create returns null).
    pub fn new(config: &Config) -> Result<Option<Self>> {
        let args = &config.args;
        let Some(msi_file) = args.value("ref_genome_msi_file") else { return Ok(None) };
        // read and checked as Java does; REDUX 2.1 uses none of these three values (R13). Corrected: the minimum
        // mapping quality is applied to the reads of the jitter analysis; the other two stay without effect
        let min_map_quality = args.integer("msi_min_map_quality", DEFAULT_MIN_MAPPING_QUALITY)?;
        args.integer("msi_max_sites_per_type", DEFAULT_NUM_SITES_PER_TYPE)?;
        args.decimal("max_site_alt_contribution", DEFAULT_MAX_SINGLE_SITE_ALT_CONTRIBUTION)?;
        // ConfigBuilder.addPath: a path that is given must exist
        for name in ["ref_genome_msi_file", "msi_model_coefficients", "msi_model_error_rates"] {
            if let Some(path) = args.value(name) {
                if !std::path::Path::new(path).exists() {
                    bail!("invalid path for -{name}: {path}");
                }
            }
        }
        let jitter_config = JitterConfig {
            sample_id: config.sample_id.clone(),
            output_dir: config.output_dir.clone(),
            sequencing_type: config.sequencing_type,
            uses_duplex_umis: config.umi_enabled && config.umi_duplex,
            apply_bq_filter: !args.has_flag("msi_no_bq_filter"),
            min_map_quality: if crate::compat::java("R13") { DEFAULT_MIN_MAPPING_QUALITY } else { min_map_quality },
            write_site_file: args.has_flag("write_msi_site_file"),
            coefficients_file: args.value("msi_model_coefficients").map(str::to_string),
            error_rates_file: args.value("msi_model_error_rates").map(str::to_string),
        };
        if args.has_flag("write_msi_plots") {
            crate::warn!("write_msi_plots is not supported: no MSI jitter plot will be written");
        }
        let sites = Sites::load(msi_file, config)?;
        crate::info!("loaded {} microsatellites regions", sites.len());
        Ok(Some(Self::from_parts(jitter_config, sites)))
    }

    fn from_parts(config: JitterConfig, sites: Sites) -> Self {
        let totals = Mutex::new(Totals::new(sites.len()));
        MsJitterAnalyser { config, sites, totals }
    }

    /// A handle for one worker thread.
    pub fn processor(&self) -> MsJitterProcessor<'_> {
        MsJitterProcessor {
            shared: self,
            min_map_quality: self.config.min_map_quality.clamp(0, 255) as u8,
            is_sbx: self.config.sequencing_type == SequencingType::Sbx,
            apply_bq_filter: self.config.apply_bq_filter,
            store_lengths: self.config.write_site_file,
            chromosome: String::new(),
            range: None,
            cursor: 0,
            window_base: 0,
            window: Vec::new(),
            types_used: [false; 3],
            far_jitter: HashMap::new(),
            repeat_lengths: HashMap::new(),
        }
    }

    /// MsJitterAnalyser.writeAnalysisOutput, called once after all reads are processed and all processors are dropped.
    pub fn write_analysis_output(&self) -> Result<()> {
        let totals = self.totals.lock().unwrap();
        let consensus_types = self.config.consensus_types();

        // MicrosatelliteSiteData.shouldKeepSite
        let java_total = crate::compat::java("R11");
        let real_variant: Vec<bool> = (0..self.sites.len()).map(|i| totals.is_real_variant_with(i, java_total)).collect();

        if self.config.write_site_file {
            let filename = format!("{}/{}.repeat.tsv.gz", self.config.output_dir, self.config.sample_id);
            let mut out = table::create_writer(&filename)?;
            self.write_site_file(&mut out, &totals, &real_variant, &consensus_types)?;
            out.flush()?;
        }

        // the nine chart tables of each consensus type
        let mut chart_tables: Vec<JitterCountsTable> = Vec::new();
        for &ct in &consensus_types {
            chart_tables.extend(self.build_group_tables(&totals, &real_variant, ct));
        }
        let filename = table::counts_table_filename(&self.config.output_dir, &self.config.sample_id);
        let mut out = table::create_writer(&filename)?;
        table::write_counts_tables(&mut out, &chart_tables)?;
        out.flush()?;
        drop(out);

        if let Some((coefficients_file, error_rates_file)) = self.config.msi_prediction_files() {
            let params = ms_model::MsModelParams::for_sequencing_type(self.config.sequencing_type);
            let calculator = ms_model::MsModelCalculator::from_files(params, coefficients_file, error_rates_file)?;
            let predicted = calculator.calc_predicted_msi_rate(&chart_tables)?;
            table::write_msi_prediction(&table::msi_prediction_filename(&self.config.output_dir, &self.config.sample_id), predicted)?;
        }

        let fitted = fit_jitter_models(&chart_tables)?;
        let filename = table::model_params_filename(&self.config.output_dir, &self.config.sample_id);
        let mut out = table::create_writer(&filename)?;
        table::write_model_params(&mut out, &fitted)?;
        out.flush()?;
        Ok(())
    }

    /// MsJitterAnalyser.buildJitterCountsTable for the nine groups of one consensus type, in one pass over the sites.
    /// Java builds each table twice looking for outliers, but never marks one, so the second pass repeats the first.
    fn build_group_tables(&self, totals: &Totals, real_variant: &[bool], ct: ConsensusType) -> Vec<JitterCountsTable> {
        let mut tables: Vec<JitterCountsTable> = (0..GROUP_COUNT).map(|g| JitterCountsTable::new(&group_name(g), ct)).collect();
        for site in 0..self.sites.len() {
            let group = self.sites.group[site];
            if group == NO_GROUP || real_variant[site] {
                continue;
            }
            // a kept site adds its row even when it has no reads of this consensus type
            let row = tables[group as usize].get_or_create_row(self.sites.repeat_count[site]);
            totals.for_each_jitter(site, ct, |jitter, count| row.add_reads(jitter, count));
        }
        tables
    }

    /// MicrosatelliteSiteFile.write: every site for each consensus type, sorted by chromosome name, start and end.
    fn write_site_file(&self, out: &mut dyn Write, totals: &Totals, real_variant: &[bool], consensus_types: &[ConsensusType]) -> Result<()> {
        let sites = &self.sites;
        let mut chromosome_order: Vec<usize> = (0..sites.chromosomes.len()).collect();
        chromosome_order.sort_by(|a, b| sites.chromosomes[*a].0.cmp(&sites.chromosomes[*b].0));
        let mut order: Vec<usize> = Vec::with_capacity(sites.len());
        for c in chromosome_order {
            let (lo, hi) = (sites.chromosomes[c].1, sites.chromosomes[c].2);
            let from = order.len();
            order.extend(lo..hi);
            // the Java sort is stable from file order
            order[from..].sort_by_key(|&i| (sites.start[i], sites.end[i], sites.file_index[i]));
        }

        let mut columns: Vec<String> =
            ["chromosome", "start", "end", "unit", "consensusType", "numReads", "numReadsRejected", "realVariant"].iter().map(|s| s.to_string()).collect();
        columns.extend(table::count_columns());
        columns.push("readRepeatLengths".to_string());
        writeln!(out, "{}", columns.join("\t"))?;

        for &ct in consensus_types {
            for &i in &order {
                write!(
                    out,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    sites.chromosome_of(i),
                    sites.start[i],
                    sites.end[i],
                    String::from_utf8_lossy(sites.unit(i)),
                    ct.name(),
                    totals.read_count(i, ct),
                    totals.rejected_count(i, ct),
                    real_variant[i]
                )?;
                for jitter in -10..=10 {
                    write!(out, "\t{}", totals.jitter_count(i, ct, jitter))?;
                }
                // MicrosatelliteSiteData.allPassingRepeats: each repeat length once for every read, ascending
                out.write_all(b"\t")?;
                let mut first = true;
                for (&(_, _, length), &count) in totals.repeat_lengths.range((i as u32, ct as u8, i32::MIN)..=(i as u32, ct as u8, i32::MAX)) {
                    for _ in 0..count {
                        if !first {
                            out.write_all(b",")?;
                        }
                        first = false;
                        write!(out, "{length}")?;
                    }
                }
                out.write_all(b"\n")?;
            }
        }
        Ok(())
    }
}

/// MsJitterAnalyser.fitJitterModels: for each consensus type the six named groups and the 3-5 base units together.
/// `chart_tables` holds the nine tables of each type; the first six are the tables Java rebuilds here, and the
/// 3-5 base table is the sum of the last three.
fn fit_jitter_models(chart_tables: &[JitterCountsTable]) -> Result<Vec<JitterModelParams>> {
    let grid = SearchGrid::new();
    let mut fitted = Vec::new();
    for tables in chart_tables.chunks(GROUP_COUNT) {
        for t in &tables[..6] {
            fitted.push(fit::perform_fit(&grid, t)?);
        }
        let mut merged = JitterCountsTable::new("3-5bp repeat", tables[0].consensus_type);
        for t in &tables[6..] {
            merged.merge_table(t);
        }
        fitted.push(fit::perform_fit(&grid, &merged)?);
    }
    Ok(fitted)
}

/// The per-thread part of SampleReadProcessor and MicrosatelliteSiteData.addReadToStats.
pub struct MsJitterProcessor<'a> {
    shared: &'a MsJitterAnalyser,
    /// the lowest mapping quality of a read that is analysed
    pub min_map_quality: u8,
    is_sbx: bool,
    apply_bq_filter: bool,
    store_lengths: bool,
    /// chromosome of the last read and its range of site indices
    chromosome: String,
    range: Option<(usize, usize)>,
    /// the first site of the chromosome that starts at or after the last read's start
    cursor: usize,
    /// counts of the sites window_base.., SITE_STRIDE values each
    window_base: usize,
    window: Vec<u32>,
    types_used: [bool; 3],
    far_jitter: HashMap<(u32, u8, i32), u32>,
    repeat_lengths: HashMap<(u32, u8, i32), u32>,
}

impl MsJitterProcessor<'_> {
    /// MsJitterAnalyser.processRead. `chromosome` is the record's reference name.
    #[inline]
    pub fn process_read(&mut self, rec: &bam::Record, chromosome: &str) {
        if rec.is_unmapped() {
            return;
        }
        if self.chromosome != chromosome {
            self.set_chromosome(chromosome);
        }
        let Some((lo, hi)) = self.range else { return };
        let sites = &self.shared.sites;
        let starts = &sites.start;
        let read_start = rec.pos() as i32 + 1;

        // move the cursor to the first site starting at or after the read start: reads arrive sorted, so this is
        // usually zero or one step, with a binary search after a jump
        let mut c = self.cursor;
        if c > lo && starts[c - 1] >= read_start {
            c = lo + starts[lo..c].partition_point(|&s| s < read_start);
        } else {
            let mut steps = 0;
            while c < hi && starts[c] < read_start {
                c += 1;
                steps += 1;
                if steps == 8 {
                    c += starts[c..hi].partition_point(|&s| s < read_start);
                    break;
                }
            }
        }
        self.cursor = c;
        if c >= hi {
            return;
        }

        let cigar = rec.raw_cigar();
        let read_end = read_start + read::reference_length(cigar) - 1;
        if starts[c] > read_end {
            return;
        }

        // ImmutableIntervalTree.containedIntervals: sites that start and end within the alignment
        let mut base_type: Option<ConsensusType> = None;
        let mut i = c;
        while i < hi && starts[i] <= read_end {
            if sites.end[i] <= read_end {
                self.add_read_to_stats(i, rec, read_start, read_end, cigar, &mut base_type);
            }
            i += 1;
        }
    }

    #[cold]
    fn set_chromosome(&mut self, chromosome: &str) {
        self.chromosome.clear();
        self.chromosome.push_str(chromosome);
        self.range = self.shared.sites.chromosome_range(chromosome).filter(|r| r.0 < r.1);
        self.cursor = self.range.map(|r| r.0).unwrap_or(0);
    }

    /// MicrosatelliteSiteData.addReadToStats
    #[inline]
    fn add_read_to_stats(&mut self, site_index: usize, rec: &bam::Record, read_start: i32, read_end: i32, cigar: &[u32], base_type: &mut Option<ConsensusType>) {
        let site = self.shared.sites.site_ref(site_index);
        let ms_read = read::analyse(&site, rec, read_start, read_end, cigar, self.apply_bq_filter);

        // the UT attribute is read once for a read, and only when a site needs it
        let base = *base_type.get_or_insert_with(|| read::extract_consensus_type(rec));
        let consensus_type = if self.is_sbx { read::sbx_consensus_type(base, &site, rec, read_start, cigar) } else { base };
        let t = consensus_type as usize;

        let offset = self.window_offset(site_index);
        let type_offset = offset + 1 + t * TYPE_STRIDE;
        if !ms_read.valid {
            self.window[type_offset + 1] += 1;
            self.types_used[t] = true;
            return;
        }
        self.window[offset] += 1;

        // SBX: a read formed from simplex reads or with medium quality bases counts as passing and nothing else
        if self.is_sbx && (consensus_type == ConsensusType::Single || ms_read.has_medium_quals) {
            return;
        }
        self.types_used[t] = true;
        self.window[type_offset] += 1;
        if (-DENSE_JITTER..=DENSE_JITTER).contains(&ms_read.jitter) {
            self.window[type_offset + 2 + (ms_read.jitter + DENSE_JITTER) as usize] += 1;
        } else {
            *self.far_jitter.entry((site_index as u32, t as u8, ms_read.jitter)).or_insert(0) += 1;
        }
        if self.store_lengths {
            *self.repeat_lengths.entry((site_index as u32, t as u8, ms_read.repeat_length)).or_insert(0) += 1;
        }
    }

    /// The offset of a site's counts in the window, moving or growing the window as needed.
    #[inline]
    fn window_offset(&mut self, site_index: usize) -> usize {
        if self.window.is_empty() {
            self.window_base = site_index;
        } else if site_index < self.window_base || site_index - self.window_base >= MAX_WINDOW_SITES {
            self.flush();
            self.window_base = site_index;
        }
        let offset = (site_index - self.window_base) * SITE_STRIDE;
        if offset + SITE_STRIDE > self.window.len() {
            self.window.resize(offset + SITE_STRIDE, 0);
        }
        offset
    }

    /// Adds this worker's counts to the shared totals.
    fn flush(&mut self) {
        if self.window.is_empty() && self.far_jitter.is_empty() && self.repeat_lengths.is_empty() {
            return;
        }
        let mut totals = self.shared.totals.lock().unwrap();
        let site_count = totals.passing.len();
        let window_sites = self.window.len() / SITE_STRIDE;
        for k in 0..window_sites {
            let passing = self.window[k * SITE_STRIDE];
            if passing != 0 {
                let p = &mut totals.passing[self.window_base + k];
                *p = p.wrapping_add(passing);
            }
        }
        for t in 0..3 {
            if !self.types_used[t] {
                continue;
            }
            let values = &mut totals.by_type[t];
            if values.is_empty() {
                values.resize(site_count * TYPE_STRIDE, 0);
            }
            let target = &mut values[self.window_base * TYPE_STRIDE..(self.window_base + window_sites) * TYPE_STRIDE];
            for k in 0..window_sites {
                let from = &self.window[k * SITE_STRIDE + 1 + t * TYPE_STRIDE..][..TYPE_STRIDE];
                for (a, b) in target[k * TYPE_STRIDE..][..TYPE_STRIDE].iter_mut().zip(from) {
                    *a = a.wrapping_add(*b);
                }
            }
        }
        for (key, count) in self.far_jitter.drain() {
            let c = totals.far_jitter.entry(key).or_insert(0);
            *c = c.wrapping_add(count);
        }
        for (key, count) in self.repeat_lengths.drain() {
            let c = totals.repeat_lengths.entry(key).or_insert(0);
            *c = c.wrapping_add(count);
        }
        drop(totals);
        self.window.clear();
        self.types_used = [false; 3];
    }
}

impl Drop for MsJitterProcessor<'_> {
    fn drop(&mut self) {
        self.flush();
    }
}

#[cfg(test)]
mod tests;
