//! ESVEE caller step (Java package esvee.caller, entry point CallerApplication).
//!
//! Reads the breakends of the depth-annotated VCF, pairs them into variants, applies the filters,
//! the panel of normals and the repeat mask, and writes the unfiltered, somatic and germline VCFs.

mod annotate;
mod config;
#[cfg(test)]
pub(crate) use config::accepted_options;
mod filters;
mod pon;
#[cfg(test)]
mod tests;
mod types;
mod variant;
mod vcf;

use crate::common::compat;
use crate::{debug, error, info, warn};
use annotate::{HotspotCache, RepeatMask, TargetRegions};
use anyhow::{bail, Context, Result};
use config::*;
use filters::*;
use pon::PonCache;
use std::io::Write;
use types::*;
use variant::*;
use vcf::*;

const HOTSPOT: &str = "HOTSPOT";
const PON_COUNT: &str = "PON_COUNT";
const REPEAT_MASK_REPEAT_CLASS: &str = "INSRMRC";
const REPEAT_MASK_REPEAT_TYPE: &str = "INSRMRT";
const REPEAT_MASK_COVERAGE: &str = "INSRMP";

const FRAG_LENGTH_DIST_PERCENTILE: f64 = 0.9975;
const FRAG_LENGTH_1_STD_DEV_PERCENTILE: f64 = 0.16;

const SV_GERMLINE_AF_THRESHOLD: f64 = 0.1;
const SV_GERMLINE_AD_THRESHOLD: f64 = 0.01;

/// log the message and stop, as Java does with System.exit(1)
fn exit_with(message: &str) -> ! {
    error!("{message}");
    std::process::exit(1);
}

/// FragmentSizeDistribution.loadFragmentLengthBounds
fn load_fragment_length_bounds(filename: &str) -> Result<FragmentLengthBounds> {
    let text = std::fs::read_to_string(filename).with_context(|| format!("failed to read fragment length file: {filename}"))?;
    let mut frequencies: Vec<(i32, i32)> = Vec::new();
    for line in text.lines().skip(1) {
        let mut values = line.split('\t');
        let parse = |v: Option<&str>| v.and_then(java_parse_int).with_context(|| format!("invalid fragment length line: {line}"));
        let length = parse(values.next())?;
        let frequency = parse(values.next())?;
        frequencies.push((length, frequency));
    }
    Ok(calculate_fragment_length_bounds(&frequencies))
}

/// FragmentSizeDistribution.calculateFragmentLengthBounds
pub(crate) fn calculate_fragment_length_bounds(frequencies: &[(i32, i32)]) -> FragmentLengthBounds {
    if frequencies.is_empty() {
        return FragmentLengthBounds { lower_bound: 0, upper_bound: 0, median: 0, std_deviation: 0.0 };
    }
    let total_fragments: i32 = frequencies.iter().fold(0i32, |a, f| a.wrapping_add(f.1));
    let mut cumulative_total = 0;
    let required_min_total = (total_fragments as f64 * (1.0 - FRAG_LENGTH_DIST_PERCENTILE)).floor() as i32;
    let required_std_dev_total = (total_fragments as f64 * FRAG_LENGTH_1_STD_DEV_PERCENTILE).floor() as i32;
    let required_max_total = (total_fragments as f64 * FRAG_LENGTH_DIST_PERCENTILE).floor() as i32;
    let median_fragment = total_fragments / 2;

    let (mut lower_bound, mut upper_bound, mut std_dev_length, mut median) = (0, 0, 0, 0);
    for &(length, frequency) in frequencies {
        if lower_bound == 0 && frequency + cumulative_total >= required_min_total {
            lower_bound = length;
        }
        if std_dev_length == 0 && frequency + cumulative_total >= required_std_dev_total {
            std_dev_length = length;
        }
        if median == 0 && frequency + cumulative_total >= median_fragment {
            median = length;
        }
        if upper_bound == 0 && frequency + cumulative_total >= required_max_total {
            upper_bound = length;
        } else {
            cumulative_total += frequency;
        }
    }
    FragmentLengthBounds { lower_bound, upper_bound, median, std_deviation: (median - std_dev_length) as f64 }
}

/// EsveeDiscordantStats.read: the second line holds the total reads, the prep reads and a count per
/// discordant fragment type
fn load_discordant_stats(filename: &str) -> Result<DiscordantStats> {
    let text = std::fs::read_to_string(filename)?;
    let line = text.lines().nth(1).context("discordant stats file without a data line")?;
    let values: Vec<&str> = line.split('\t').collect();
    let parse = |i: usize| -> Result<i64> {
        let v = values.get(i).with_context(|| format!("discordant stats line with {} fields", values.len()))?;
        match v.parse::<i64>() {
            Ok(value) => Ok(value),
            Err(e) if compat::java("C12") => Err(e).with_context(|| format!("NumberFormatException: For input string: \"{v}\"")),
            Err(_) => bail!("line 2, field {}: '{v}' is not a whole number", i + 1),
        }
    };
    let mut stats = DiscordantStats { total_reads: parse(0)?, ..Default::default() };
    parse(1)?; // the reads written by prep, which the caller does not use
    for i in 0..DISCORDANT_FRAG_TYPES {
        stats.type_counts[i] = parse(i + 2)?;
    }
    Ok(stats)
}

/// CallerApplication.isGermline
pub(crate) fn is_germline(var: &Variant, reference_id: Option<&str>) -> bool {
    // a match to a known germline variant
    if var.is_saga_matched {
        return true;
    }
    let breakend = var.start();
    let (mut max_germline_af, mut max_tumor_af) = (0.0f64, 0.0f64);
    let (mut germline_ad, mut tumor_ad) = (0, 0);
    for genotype in &breakend.context.genotypes {
        let af = var.calc_allelic_frequency(genotype);
        if reference_id == Some(genotype.sample.as_str()) {
            max_germline_af = max_germline_af.max(af);
            germline_ad = Breakend::fragment_count_of(genotype);
        } else {
            max_tumor_af = max_tumor_af.max(af);
            tumor_ad = Breakend::fragment_count_of(genotype);
        }
    }
    if max_germline_af >= SV_GERMLINE_AF_THRESHOLD * max_tumor_af {
        // also check the relative fragment counts
        let ad_ratio = if tumor_ad > 0 { germline_ad as f64 / tumor_ad as f64 } else { 1.0 };
        if ad_ratio >= SV_GERMLINE_AD_THRESHOLD {
            return true;
        }
    }
    false
}

/// CallerApplication.processVariant, the part for -manual_ref_depth: a sample without a depth gets
/// the given one, for the reads and for the read pairs.
pub(crate) fn apply_manual_ref_depth(record: &mut VcfRecord, depth: i32, genotype_ids: &GenotypeIds, java_c5: bool) -> Result<()> {
    let set = |genotype: &mut Genotype| {
        if !genotype.has_attr(REF_DEPTH) {
            genotype.set_attr(REF_DEPTH, depth.to_string());
            genotype.set_attr(REF_DEPTH_PAIR, depth.to_string());
        }
    };
    if java_c5 {
        // Java (C5): only done when there is a reference sample, and then for the tumor as well, which
        // leaves a tumor-only run without values and fails in a reference-only run
        if genotype_ids.has_reference() {
            for sample in [&genotype_ids.reference_id, &genotype_ids.tumor_id] {
                let Some(genotype) = record.genotype_mut(sample) else {
                    bail!("NullPointerException: manual ref depth without a genotype for sample '{sample}'");
                };
                set(genotype);
            }
        }
        return Ok(());
    }
    // each sample of the run
    for (in_run, sample) in [(genotype_ids.has_reference(), &genotype_ids.reference_id), (genotype_ids.has_tumor(), &genotype_ids.tumor_id)] {
        if !in_run {
            continue;
        }
        let Some(genotype) = record.genotype_mut(sample) else {
            bail!("manual ref depth: record {} at {}:{} has no genotype for sample '{sample}'", record.id, record.chrom, record.pos);
        };
        set(genotype);
    }
    Ok(())
}

/// SvDataCache: pairs breakend records into variants as they are read
struct SvDataCache<'a> {
    config: &'a CallerConfig,
    target_regions: &'a TargetRegions,
    reference_ordinal: i32,
    tumor_ordinal: i32,
    variants: Vec<Variant>,
    unmatched: hashbrown::HashMap<String, VcfRecord>,
    java_c14: bool,
    /// breakends on a contig other than 1-22, X and Y, and breakends left out because their mate is
    /// on such a contig (the latter only counted when C14 is corrected)
    other_contig: usize,
    mate_on_other_contig: usize,
    /// breakends on chromosomes 1-22, X and Y, and how many of them are named as the genome version
    /// names them; the first name that is not
    human_contig: usize,
    version_named: usize,
    first_misnamed: Option<String>,
}

/// The contig of the mate of a paired breakend, from its ALT (`A[chr2:100[`)
pub(crate) fn mate_contig(alt: &str) -> Option<&str> {
    let (_, _, mate, _) = match_breakend(alt)?;
    mate.rsplit_once(':').map(|(contig, _)| contig)
}

impl SvDataCache<'_> {
    fn new<'a>(config: &'a CallerConfig, target_regions: &'a TargetRegions, genotype_ids: &GenotypeIds) -> SvDataCache<'a> {
        SvDataCache {
            config,
            target_regions,
            reference_ordinal: genotype_ids.reference_ordinal,
            tumor_ordinal: genotype_ids.tumor_ordinal,
            variants: Vec::new(),
            unmatched: hashbrown::HashMap::new(),
            java_c14: compat::java("C14"),
            other_contig: 0,
            mate_on_other_contig: 0,
            human_contig: 0,
            version_named: 0,
            first_misnamed: None,
        }
    }

    fn process_variant(&mut self, record: VcfRecord) {
        if !is_human_chromosome(&record.chrom) {
            self.other_contig += 1;
            return;
        }
        self.human_contig += 1;
        if record.chrom.starts_with(CHR_PREFIX) != self.config.ref_gen_version.is37() {
            self.version_named += 1;
        } else if self.first_misnamed.is_none() {
            self.first_misnamed = Some(record.chrom.clone());
        }
        if is_single_breakend(&record) {
            let mut variant = Variant::new_sgl(record, self.reference_ordinal, self.tumor_ordinal);
            if self.target_regions.has_target_regions()
                && !self.target_regions.in_target_regions(variant.chromosome_start(), variant.pos_start(), PANEL_INCLUSION_BUFFER)
            {
                variant.add_filter(FilterType::OffTarget);
            }
            self.add_sv_data(variant);
            return;
        }
        let Some(mate_id) = record.attr_string_opt(MATE_ID) else { return };

        // Java (C14): the mate on another contig was skipped above, so this breakend waits for it for
        // ever and is only seen in the count of unmatched breakends. Corrected: the variant is left
        // out as a whole, knowingly, and reported as such (see the notes: output for chromosomes
        // 1-22, X and Y only, and no breakend without its mate).
        if !self.java_c14 && mate_contig(record.alt()).is_some_and(|contig| !contig.is_empty() && !is_human_chromosome(contig)) {
            self.mate_on_other_contig += 1;
            return;
        }

        // the variant is complete once both breakends have been read
        let Some(first) = self.unmatched.remove(&mate_id) else {
            self.unmatched.insert(record.id.clone(), record);
            return;
        };
        let mut variant = Variant::new_sv(first, record, self.reference_ordinal, self.tumor_ordinal);

        // one of the breakends at least must be within a targeted region
        if self.target_regions.has_target_regions() {
            let within = |b: &Breakend| self.target_regions.in_target_regions(&b.chromosome, b.position, PANEL_INCLUSION_BUFFER);
            if !within(variant.start()) && !within(variant.breakend(SE_END)) {
                variant.add_filter(FilterType::OffTarget);
            }
        }
        self.add_sv_data(variant);
    }

    fn add_sv_data(&mut self, var: Variant) {
        // optionally keep only the variants with a breakend in one of the specified regions
        let regions = &self.config.specific_regions;
        if !regions.is_empty() {
            let within = |chromosome: &str, position: i32| regions.iter().any(|r| r.contains_position(chromosome, position));
            if !within(var.chromosome_start(), var.pos_start()) && (var.is_sgl() || !within(var.chromosome_end(), var.pos_end())) {
                return;
            }
        }
        // optionally keep only the specified chromosomes
        let chromosomes = &self.config.specific_chromosomes;
        if !chromosomes.is_empty()
            && !chromosomes.iter().any(|c| c == var.chromosome_start())
            && (var.sv_type == SvType::Sgl || !chromosomes.iter().any(|c| c == var.chromosome_end()))
        {
            return;
        }
        self.variants.push(var);
    }
}

/// The output files of a run: unfiltered always, somatic with a tumor, germline with a reference
struct OutputVcfs {
    header: VcfHeader,
    counts: HeaderCounts,
    unfiltered: VcfWriter,
    somatic: Option<VcfWriter>,
    germline: Option<VcfWriter>,
}

impl OutputVcfs {
    /// VcfWriter: the input header with the sample columns of the run and the caller's own lines
    fn create(config: &CallerConfig, input_header: &VcfHeader) -> Result<OutputVcfs> {
        let mut header = VcfHeader { meta_lines: input_header.meta_lines.clone(), samples: Vec::new() };
        if let Some(reference_id) = &config.reference_id {
            header.samples.push(reference_id.clone());
        }
        if let Some(tumor_id) = &config.tumor_id {
            header.samples.push(tumor_id.clone());
        }
        header.add_meta_line(format!("##esveeVersion={}", crate::VERSION));
        header.add_meta_line(filter_line("PASS", "Variant passes all filters"));
        for filter in FILTER_TYPES {
            header.add_meta_line(filter_line(filter.vcf_tag(), filter.vcf_desc()));
        }
        header.add_meta_line(info_line(HOTSPOT, "0", "Flag", "Variant is a hotspot"));
        header.add_meta_line(info_line(PON_COUNT, "1", "Integer", "PON count if in PON"));
        header.add_meta_line(info_line(REPEAT_MASK_REPEAT_CLASS, "1", "String", "Inserted sequence repeatmasker repeat class"));
        header.add_meta_line(info_line(REPEAT_MASK_REPEAT_TYPE, "1", "String", "Inserted sequence repeatmasker repeat type"));
        header.add_meta_line(info_line(
            REPEAT_MASK_COVERAGE,
            "1",
            "Float",
            "Portion of inserted sequence whose alignment overlaps the repeatmasker repeat",
        ));

        let sample_id = config.file_sample_id();
        let filename = |file_id: &str| {
            // SAMPLE_ID.esvee.run_id.somatic.vcf.gz, or without the output id if there is none
            let output_id = match &config.output_id {
                Some(id) => format!("{id}.{file_id}"),
                None => file_id.to_string(),
            };
            form_output_file(&config.output_dir, &sample_id, "esvee", "vcf.gz", Some(&output_id))
        };
        if header.samples.len() == 2 && header.samples[0] == header.samples[1] {
            // htsjdk rejects the header once the first file has been created, which is left empty
            std::fs::File::create(filename("unfiltered"))?;
            bail!("Your input file has a malformed header: BUG: VCF header has duplicate sample names");
        }
        let unfiltered = VcfWriter::create(&filename("unfiltered"), &header)?;
        let somatic = if config.has_tumor() { Some(VcfWriter::create(&filename("somatic"), &header)?) } else { None };
        let germline = if config.has_reference() { Some(VcfWriter::create(&filename("germline"), &header)?) } else { None };
        let counts = header.counts();
        Ok(OutputVcfs { header, counts, unfiltered, somatic, germline })
    }

    /// VcfWriter.buildGenotype: the sample's genotype with AD and DP set from the depth fields
    fn build_genotype(breakend: &Breakend, sample_id: &str) -> Result<Genotype> {
        let Some(existing) = breakend.context.genotype(sample_id) else {
            bail!("NullPointerException: sample {sample_id} has no genotype in the VCF under that name");
        };
        let mut genotype = existing.clone();
        let ref_pair_support = existing.attr_int(REF_DEPTH_PAIR, 0);
        let ref_support = existing.attr_int(REF_DEPTH, 0);
        let var_support = existing.attr_int(TOTAL_FRAGS, 0);
        let total_depth = ref_pair_support.wrapping_add(ref_support).wrapping_add(var_support);
        genotype.ad = Some(format!("{},{}", ref_pair_support.wrapping_add(ref_support), var_support));
        genotype.dp = Some(total_depth.to_string());
        Ok(genotype)
    }

    fn write_breakend(&mut self, config: &CallerConfig, genotype_ids: &GenotypeIds, var: &Variant, breakend: &Breakend) -> Result<()> {
        let mut genotypes: Vec<Genotype> = Vec::with_capacity(2);
        if genotype_ids.has_tumor() {
            genotypes.push(Self::build_genotype(breakend, config.tumor_id.as_deref().unwrap_or(""))?);
        }
        if genotype_ids.has_reference() {
            genotypes.push(Self::build_genotype(breakend, config.reference_id.as_deref().unwrap_or(""))?);
        }

        let mut record = VcfRecord {
            chrom: breakend.context.chrom.clone(),
            pos: breakend.context.pos,
            id: breakend.context.id.clone(),
            ref_allele: breakend.context.ref_allele.clone(),
            alts: breakend.context.alts.clone(),
            log10_p_error: breakend.context.phred_scaled_qual() / -10.0,
            filters: Vec::new(),
            info: breakend.context.info.clone(),
            genotypes: Vec::new(),
        };
        if var.is_hotspot {
            record.set_attr(HOTSPOT, None);
        }
        if var.pon_count > 0 {
            record.set_attr(PON_COUNT, Some(var.pon_count.to_string()));
        }
        if let Some(rm) = &var.rm_annotation {
            record.set_attr(REPEAT_MASK_REPEAT_CLASS, Some(rm.class_type.clone()));
            record.set_attr(REPEAT_MASK_REPEAT_TYPE, Some(rm.repeat.clone()));
            record.set_attr(REPEAT_MASK_COVERAGE, Some(crate::common::fmt::vcf_double(rm.coverage)));
        }
        let encoded = encode_record(&self.header, &self.counts, &record, record.log10_p_error, &genotypes)?;

        // first the unfiltered VCF with all breakends
        let all_filters = var.filters;
        let mut tags: Vec<&str> = all_filters.iter().map(|f| f.vcf_tag()).collect();
        self.unfiltered.write(&encoded, &mut tags)?;

        if var.germline {
            let germline_filters: Vec<FilterType> = all_filters.iter().filter(|f| !f.germline_only()).collect();
            if let Some(writer) = &mut self.germline {
                if germline_filters.is_empty() {
                    writer.write(&encoded, &mut Vec::new())?;
                }
            }
        } else if let Some(writer) = &mut self.somatic {
            if all_filters.is_empty() || (all_filters.len() == 1 && all_filters.contains(FilterType::Pon)) {
                writer.write(&encoded, &mut tags)?;
            }
        }
        Ok(())
    }

    /// A failure while writing stops Java with its files created but still empty.
    fn abort(self) {
        let paths: Vec<String> =
            [Some(&self.unfiltered), self.somatic.as_ref(), self.germline.as_ref()].into_iter().flatten().map(|w| w.path().to_string()).collect();
        drop(self);
        for path in paths {
            let _ = std::fs::File::create(path);
        }
    }

    fn close(self) -> Result<()> {
        self.unfiltered.close()?;
        if let Some(writer) = self.somatic {
            writer.close()?;
        }
        if let Some(writer) = self.germline {
            writer.close()?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------------------------
// breakend TSV (BreakendWriter)

struct BreakendData {
    breakend_id: String,
    chromosome: String,
    position: i32,
    orient: Orient,
    sv_type: SvType,
    mate_chromosome: String,
    mate_position: i32,
    mate_orient: Orient,
    insert_sequence: String,
    values_str: String,
}

impl BreakendData {
    fn matches(&self, variants: &[Variant], other: BreakendRef) -> bool {
        let var = &variants[other.var as usize];
        let breakend = var.breakend(other.se as usize);
        if self.breakend_id != breakend.vcf_id() {
            return false;
        }
        if self.sv_type != var.sv_type {
            // allow DUP vs INS differences
            let dup_ins = (self.sv_type == SvType::Dup && var.sv_type == SvType::Ins) || (self.sv_type == SvType::Ins && var.sv_type == SvType::Dup);
            if !dup_ins {
                return false;
            }
        }
        if self.chromosome != breakend.chromosome || self.position != breakend.position || self.orient != breakend.orient {
            return false;
        }
        if self.sv_type == SvType::Sgl {
            return self.insert_sequence == breakend.insert_sequence;
        }
        let mate = var.breakend(1 - other.se as usize);
        self.mate_chromosome == mate.chromosome && self.mate_position == mate.position && self.mate_orient == mate.orient
    }
}

/// FileWriterUtils.pathFromFile
fn path_from_file(file_path: &str) -> String {
    match file_path.rfind('/') {
        Some(i) => {
            let parent = file_path[..i].trim_end_matches('/');
            if parent.is_empty() { "//".to_string() } else { format!("{parent}/") }
        }
        None => "./".to_string(),
    }
}

/// BreakendWriter.loadBreakendData: the breakends of the assembly step's TSV by chromosome, sorted
fn load_breakend_data(text: &str) -> Result<(String, Vec<(String, Vec<BreakendData>)>)> {
    let mut lines = text.lines();
    let header_line = lines.next().unwrap_or("").to_string();
    let columns: Vec<&str> = header_line.split('\t').collect();
    let column = |name: &str| columns.iter().rposition(|c| *c == name).with_context(|| format!("NullPointerException: breakend TSV has no column {name}"));
    let id_index = column("Id")?;
    let chr_index = column("Chromosome")?;
    let pos_index = column("Position")?;
    let orient_index = column("Orientation")?;
    let mate_chr_index = column("MateChr")?;
    let mate_pos_index = column("MatePos")?;
    let mate_orient_index = column("MateOrient")?;
    let type_index = column("Type")?;
    let ins_seq_index = column("InsertedBases")?;

    let mut map: Vec<(String, Vec<BreakendData>)> = Vec::new();
    let mut current: Option<usize> = None;
    for line in lines {
        let values: Vec<&str> = line.split('\t').collect();
        let field = |i: usize| values.get(i).copied().with_context(|| format!("ArrayIndexOutOfBoundsException: breakend TSV line with {} fields", values.len()));
        let int = |s: &str| java_parse_int(s).with_context(|| format!("NumberFormatException: For input string: \"{s}\""));
        let chromosome = field(chr_index)?;
        let sv_type = SvType::value_of(field(type_index)?);
        let (mut mate_chromosome, mut mate_position, mut mate_orient) = (String::new(), -1, Orient::Forward);
        if sv_type != SvType::Sgl {
            mate_chromosome = field(mate_chr_index)?.to_string();
            mate_position = int(field(mate_pos_index)?)?;
            mate_orient = Orient::from_byte_str(field(mate_orient_index)?);
        }
        let data = BreakendData {
            breakend_id: field(id_index)?.to_string(),
            chromosome: chromosome.to_string(),
            position: int(field(pos_index)?)?,
            orient: Orient::from_byte_str(field(orient_index)?),
            sv_type,
            mate_chromosome,
            mate_position,
            mate_orient,
            insert_sequence: field(ins_seq_index)?.to_string(),
            values_str: line.to_string(),
        };
        if current.is_none_or(|i| map[i].0 != chromosome) {
            current = Some(match map.iter().position(|(c, _)| c == chromosome) {
                Some(i) => i,
                None => {
                    map.push((chromosome.to_string(), Vec::new()));
                    map.len() - 1
                }
            });
        }
        map[current.unwrap()].1.push(data);
    }
    for (_, list) in map.iter_mut() {
        // CommonUtils.compareJunctions within a chromosome
        list.sort_by(|a, b| {
            a.position.cmp(&b.position).then_with(|| {
                if a.orient == b.orient {
                    std::cmp::Ordering::Equal
                } else if a.orient.is_forward() {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Greater
                }
            })
        });
    }
    Ok((header_line, map))
}

/// BreakendWriter.writeBreakends: the assembly step's breakend TSV with the caller's columns added
fn write_breakend_tsv(config: &CallerConfig, variants: &[Variant], map: &BreakendMap) -> Result<()> {
    let input_dir = path_from_file(&config.vcf_file);
    let sample_id = config.file_sample_id();
    let mut initial_tsv = form_output_file(&input_dir, &sample_id, "esvee", "breakend.tsv", config.output_id.as_deref());
    if !std::path::Path::new(&initial_tsv).exists() && config.output_id.is_some() {
        initial_tsv = form_output_file(&input_dir, &sample_id, "esvee", "breakend.tsv", None);
    }
    if !std::path::Path::new(&initial_tsv).exists() {
        warn!("cannot write breakend TSV since input TSV({initial_tsv}) not found");
        return Ok(());
    }
    let new_tsv = form_output_file(&config.output_dir, &sample_id, "esvee", "final_breakend.tsv", config.output_id.as_deref());

    let text = match std::fs::read_to_string(&initial_tsv) {
        Ok(t) => t,
        Err(e) => {
            error!("failed to rewrite breakends: {e}");
            return Ok(());
        }
    };
    let (header_line, data_map) = load_breakend_data(&text)?;
    if data_map.is_empty() {
        return Ok(());
    }
    debug!("annotating {} existing breakends", data_map.iter().map(|(_, l)| l.len()).sum::<usize>());

    let mut out = std::io::BufWriter::new(std::fs::File::create(&new_tsv).with_context(|| format!("creating {new_tsv}"))?);
    writeln!(out, "{header_line}\tGermline\tFilter\tVAF\tPONCount")?;

    for chromosome in human_chromosomes() {
        let chr_str = config.ref_gen_version.versioned_chromosome(&chromosome);
        let (Some(breakend_list), Some((_, data_list))) = (map.get(&chr_str), data_map.iter().find(|(c, _)| *c == chr_str)) else { continue };
        let mut data_index: usize = 0;
        for &b in breakend_list {
            if data_index >= data_list.len() {
                break;
            }
            let position = variants[b.var as usize].breakend(b.se as usize).position;
            let mut data = &data_list[data_index];
            if !data.matches(variants, b) {
                // search backwards then forwards among the entries at this position
                let mut match_found = false;
                for search_back in [true, false] {
                    let mut test_index = if search_back { data_index as i64 - 1 } else { data_index as i64 + 1 };
                    while test_index >= 0 && (test_index as usize) < data_list.len() {
                        data = &data_list[test_index as usize];
                        if data.position > position && !search_back {
                            break;
                        } else if data.position < position && search_back {
                            break;
                        }
                        if data.matches(variants, b) {
                            data_index = test_index as usize;
                            match_found = true;
                            break;
                        }
                        test_index += if search_back { -1 } else { 1 };
                    }
                    if match_found {
                        break;
                    }
                }
            }
            if data.matches(variants, b) {
                let var = &variants[b.var as usize];
                let filters_str =
                    if var.filters.is_empty() { "PASS".to_string() } else { var.filters.iter().map(|f| f.vcf_tag()).collect::<Vec<_>>().join(";") };
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}",
                    data.values_str,
                    var.germline,
                    filters_str,
                    crate::common::fmt::java_fixed(var.breakend_allelic_frequency(b.se as usize), 3),
                    var.pon_count
                )?;
                data_index += 1;
            } else {
                debug!("unmatched breakend({})", describe_breakend(variants, b));
            }
        }
    }
    out.flush()?;
    Ok(())
}

/// VariantFilters.logFilterTypeCounts
fn log_filter_type_counts(variants: &[Variant]) {
    let passing = variants.iter().filter(|v| v.filters.is_empty()).count();
    info!("variants passing({}) filtered({})", passing, variants.len() - passing);
    for filter in FILTER_TYPES {
        let count = variants.iter().filter(|v| v.filters.contains(filter)).count();
        if count > 0 {
            debug!("variant filter {}: count({})", filter.java_name(), count);
        }
    }
}

fn log_messages(messages: &[(bool, String)]) {
    for (is_error, message) in messages {
        if *is_error {
            error!("{message}");
        } else if message.contains("not ordered") {
            warn!("{message}");
        } else {
            info!("{message}");
        }
    }
}

pub fn run(argv: &[String]) -> Result<()> {
    let start = std::time::Instant::now();

    let values = match ConfigValues::parse(argv) {
        Ok(v) => v,
        Err(errors) => {
            for e in errors {
                error!("{e}");
            }
            std::process::exit(1);
        }
    };
    if values.has_flag("log_debug") || matches!(values.value("log_level"), Some("DEBUG") | Some("TRACE") | Some("ALL")) {
        crate::common::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
    }
    if crate::version_log_due() {
        info!("Esvee version {}", crate::VERSION);
    }

    let config = CallerConfig::new(&values)?;
    let filter_constants = filter_constants(&values, config.ref_gen_version);
    let version = config.ref_gen_version;

    info!("loading reference data");

    // the three sets of reference data are read at the same time; their messages keep Java's order
    let (pon_margin, pon_sgl_margin) = (values.integer("pon_margin"), values.integer("pon_sgl_margin"));
    let owned = |name: &str| values.value(name).map(str::to_string);
    let (pon_sv_file, pon_sgl_file) = (owned("pon_sv_file"), owned("pon_sgl_file"));
    let (artefact_sv_file, artefact_sgl_file) = (owned("artefact_pon_sv_file"), owned("artefact_pon_sgl_file"));
    let has_artefact_pon = artefact_sv_file.is_some() || artefact_sgl_file.is_some();
    let repeat_mask_file = owned("repeat_mask_file");

    let pon_thread = std::thread::spawn(move || PonCache::new(pon_margin, pon_sgl_margin, pon_sv_file.as_deref(), pon_sgl_file.as_deref()));
    let artefact_thread = has_artefact_pon
        .then(|| std::thread::spawn(move || PonCache::new(pon_margin, pon_sgl_margin, artefact_sv_file.as_deref(), artefact_sgl_file.as_deref())));
    let java_c11 = compat::java("C11");
    let repeat_mask_thread = repeat_mask_file.clone().map(|file| {
        std::thread::spawn(move || {
            let mut repeat_mask = RepeatMask::default();
            let result = repeat_mask.load(&file, version, java_c11);
            (repeat_mask, result)
        })
    });

    let mut pon_cache = pon_thread.join().expect("PON loading thread")?;
    pon_cache.set_specific_sgl_regions(version);
    log_messages(&pon_cache.messages);
    if !pon_cache.has_valid_data() {
        exit_with("invalid PON, exiting");
    }

    let mut artefact_pon_cache = match artefact_thread {
        Some(thread) => {
            let cache = thread.join().expect("artefact PON loading thread")?;
            log_messages(&cache.messages);
            if !cache.has_valid_data() {
                exit_with("invalid artefact PON, exiting");
            }
            Some(cache)
        }
        None => None,
    };

    let mut hotspot_cache = HotspotCache::default();
    if let Some(file) = values.value("known_hotspot_file") {
        match hotspot_cache.load(file) {
            Ok(count) => info!("loaded {count} known hotspot records from file({file})"),
            // Java (C13): a file that cannot be read is logged and the run goes on without hotspots, and
            // a bad line is an uncaught exception. Corrected: both stop the run with one message.
            Err(e) if !compat::java("C13") => bail!("failed to load hotspot data file({file}): {e:#}"),
            Err(e) => match e.downcast_ref::<std::io::Error>() {
                Some(io) => error!("failed to load hotspot data file({file}): {io}"),
                None => return Err(e),
            },
        }
    }

    let sample_id = config.file_sample_id();
    let frag_length_filename = form_prep_input_filename(&config.prep_dir, &sample_id, "fragment_length.tsv", config.output_id.as_deref());
    let disc_stats_filename = form_prep_input_filename(&config.prep_dir, &sample_id, "disc_stats.tsv", config.output_id.as_deref());
    if !std::path::Path::new(&frag_length_filename).exists() || !std::path::Path::new(&disc_stats_filename).exists() {
        if compat::java("C10") {
            // Java (C10): the two file names are passed without placeholders, so they are not shown
            exit_with("missing input files: disc-stats and frag-lengths");
        }
        let state = |file: &str| if std::path::Path::new(file).exists() { "found" } else { "missing" };
        exit_with(&format!(
            "missing input files of sample({sample_id}), which names the prep files: disc-stats({disc_stats_filename}: {}) frag-lengths({frag_length_filename}: {})",
            state(&disc_stats_filename),
            state(&frag_length_filename)
        ));
    }
    let fragment_length_bounds = load_fragment_length_bounds(&frag_length_filename)?;
    let discordant_stats = match load_discordant_stats(&disc_stats_filename) {
        Ok(stats) => stats,
        Err(e) if compat::java("C12") => {
            // Java (C12): the reader logs this and returns null, which the filters then dereference
            error!("failed to load Esvee discordant stat file({disc_stats_filename}): {e:#}");
            bail!("NullPointerException: no discordant stats");
        }
        Err(e) => bail!("failed to load Esvee discordant stat file({disc_stats_filename}): {e:#}"),
    };
    info!(
        "fragment length dist: bounds({} - {}) median({} stdDev={})",
        fragment_length_bounds.lower_bound,
        fragment_length_bounds.upper_bound,
        fragment_length_bounds.median,
        crate::common::fmt::java_fixed(fragment_length_bounds.std_deviation, 1)
    );
    let variant_filters = VariantFilters::new(filter_constants, fragment_length_bounds, discordant_stats, config.sequencing_type);

    let mut target_regions = TargetRegions::default();
    if let Some(file) = values.value("target_regions_bed") {
        let count = target_regions.load(file, version)?;
        info!("loaded {count} target regions from file({file})");
    }
    let targeted_panel_mode = target_regions.has_target_regions();

    let repeat_mask = match repeat_mask_thread {
        Some(thread) => {
            let file = repeat_mask_file.as_deref().unwrap_or("");
            let (repeat_mask, result) = thread.join().expect("repeat mask loading thread");
            match result {
                Ok(()) => info!("loaded {} repeat-mask entries from file({file})", repeat_mask.entry_count()),
                Err(e) => {
                    match e.downcast_ref::<std::io::Error>() {
                        Some(io) => error!("failed to load repeat-mask data from file({file}): {io}"),
                        None => error!("{e:#}"),
                    }
                    std::process::exit(1);
                }
            }
            repeat_mask
        }
        None => RepeatMask::default(),
    };

    // CallerApplication.run
    if config.tumor_id.is_none() && config.reference_id.is_none() {
        error!("missing sample config");
        exit_with("invalid config, exiting");
    }
    if config.vcf_file.is_empty() || !std::path::Path::new(&config.vcf_file).exists() {
        error!("missing or invalid VCF file({})", config.vcf_file);
        exit_with("invalid config, exiting");
    }

    // processVcf
    let mut input = VcfInput::open(&config.vcf_file)?;
    if config.manual_ref_depth > 0 && !input.header.has_format_line(REF_DEPTH) {
        input.header.add_meta_line(format_line(REF_DEPTH, "1", "Integer", "Count of reads mapping across this breakend"));
        input.header.add_meta_line(format_line(
            REF_DEPTH_PAIR,
            "1",
            "Integer",
            "Count of reference read pairs spanning this breakend supporting the reference allele",
        ));
    }
    info!("sample({}) processing VCF({})", sample_id, config.vcf_file);

    let genotype_ids = GenotypeIds::from_vcf_header(&input.header.samples, config.reference_id.as_deref(), config.tumor_id.as_deref());
    if !compat::java("C4") && config.has_tumor() && config.tumor_id == config.reference_id {
        // Java gets as far as the output header, which htsjdk rejects for its duplicate sample names
        exit_with(&format!("tumor and reference are the same sample({}): give two samples of the VCF, or only one of the two options", genotype_ids.tumor_id));
    }
    if (config.has_tumor() && genotype_ids.tumor_ordinal < 0) || (config.has_reference() && genotype_ids.reference_ordinal < 0) {
        if compat::java("C4") {
            exit_with(&format!("missing sample names in VCF: [{}]", input.header.samples.join(", ")));
        }
        let mut missing: Vec<String> = Vec::new();
        if config.has_tumor() && genotype_ids.tumor_ordinal < 0 {
            missing.push(format!("tumor({})", config.tumor_id.as_deref().unwrap_or("")));
        }
        if config.has_reference() && genotype_ids.reference_ordinal < 0 {
            missing.push(format!("reference({})", config.reference_id.as_deref().unwrap_or("")));
        }
        exit_with(&format!(
            "sample {} not in the VCF, whose samples are [{}]: a sample id must be the same as the VCF sample name",
            missing.join(" and "),
            input.header.samples.join(", ")
        ));
    }
    if config.germline_only() && compat::java("C9") {
        // Java (C9): the tumor's ordinal and id under ref, and the reference's under tumor
        info!(
            "germline mode ref({}: {}) tumor({}: {})",
            genotype_ids.tumor_ordinal, genotype_ids.tumor_id, genotype_ids.reference_ordinal, genotype_ids.reference_id
        );
    } else if config.germline_only() {
        info!("germline mode ref({}: {})", genotype_ids.reference_ordinal, genotype_ids.reference_id);
    } else if !config.has_reference() {
        info!("tumor genotype info({}: {})", genotype_ids.tumor_ordinal, genotype_ids.tumor_id);
    } else {
        info!(
            "genotype info: ref({}: {}) tumor({}: {})",
            genotype_ids.reference_ordinal, genotype_ids.reference_id, genotype_ids.tumor_ordinal, genotype_ids.tumor_id
        );
    }

    let mut cache = SvDataCache::new(&config, &target_regions, &genotype_ids);
    let java_c5 = compat::java("C5");
    let mut processed_variants = 0usize;
    while let Some(mut record) = input.next_record()? {
        processed_variants += 1;
        if processed_variants % 100000 == 0 {
            debug!("sample({}) processed {} variants", config.tumor_id.as_deref().unwrap_or("null"), processed_variants);
        }
        if config.manual_ref_depth > 0 {
            apply_manual_ref_depth(&mut record, config.manual_ref_depth, &genotype_ids, java_c5)?;
        }
        cache.process_variant(record);
    }
    let SvDataCache { mut variants, unmatched, other_contig, mate_on_other_contig, human_contig, version_named, first_misnamed, .. } = cache;

    // Java (C6): output goes by the chromosome names of the genome version, so breakends named the
    // other way are read, filtered and then never written; when that is all of them the three VCFs
    // come out with a header only and the run ends well. Corrected: that case stops the run.
    if !compat::java("C6") && human_contig > version_named {
        let (version_str, naming) = if version.is37() { ("37", "without") } else { ("38", "with") };
        let example = first_misnamed.as_deref().unwrap_or("");
        if version_named == 0 {
            bail!(
                "no contig of the input VCF({}) is named as ref_genome_version {version_str} names it ({naming} the chr prefix): \
                 all {human_contig} breakends are on contigs named like '{example}'; check the genome version",
                config.vcf_file
            );
        }
        warn!(
            "{} of {human_contig} breakends are on contigs not named as ref_genome_version {version_str} names them ({naming} the chr prefix), \
             such as '{example}', and are not written",
            human_contig - version_named
        );
    }
    if !compat::java("C14") && (other_contig > 0 || mate_on_other_contig > 0) {
        warn!(
            "left out {other_contig} breakends on contigs other than chromosomes 1-22, X and Y, and {mate_on_other_contig} breakends whose mate is on such a contig"
        );
    }

    if unmatched.is_empty() {
        info!(
            "loaded {} breakends, SVs({}) SGLs({})",
            processed_variants,
            variants.iter().filter(|v| !v.is_sgl()).count(),
            variants.iter().filter(|v| v.is_sgl()).count()
        );
    } else {
        warn!("loaded {} breakends with unmatched({}) complete({}) hardFiltered({})", processed_variants, unmatched.len(), variants.len(), 0);
    }
    drop(unmatched);

    info!("writing output VCF files to {}", config.output_dir);
    let mut output = OutputVcfs::create(&config, &input.header)?;

    if variants.is_empty() {
        info!("writing empty VCF");
        output.close()?;
        info!("Esvee caller complete, mins({})", crate::common::fmt::java_fixed(start.elapsed().as_millis() as f64 / 60000.0, 3));
        return Ok(());
    }

    let map = BreakendMap::build(&variants);
    mark_line_sites(&mut variants, &map);

    info!("applying filters");
    for v in 0..variants.len() {
        if hotspot_cache.is_hotspot_variant(&variants[v]) {
            variants[v].is_hotspot = true;
            // called regardless of proximity to panel boundaries
            variants[v].filters.remove(FilterType::OffTarget);
        }
        variant_filters.apply_filters(&mut variants, v);
    }
    if targeted_panel_mode {
        variant_filters.apply_adjacent_filters(&mut variants, &map);
    }

    // set germline status and final filters based on LINE
    let reference_id = config.reference_id.as_deref();
    for var in variants.iter_mut() {
        if config.germline_only() || is_germline(var, reference_id) {
            var.germline = true;
        }
    }
    for v in 0..variants.len() {
        adjust_line_sites(&mut variants, v);
    }

    deduplicate_variants(&mut variants, &map);

    if pon_cache.has_valid_data() {
        info!("applying PON filters");
        pon_cache.annotate_variants(&mut variants, &map);
    }
    if let Some(cache) = artefact_pon_cache.as_mut() {
        if cache.has_valid_data() {
            info!("applying artefacts PON filters");
            cache.annotate_variants(&mut variants, &map);
        }
    }
    if repeat_mask.has_data() {
        info!("annotating with repeat-mask information");
        let annotated = repeat_mask.annotate_variants(&mut variants);
        debug!("marked {annotated} repeat mask annotations");
    }

    // VcfWriter.writeBreakends
    let mut failure = None;
    'chromosomes: for chromosome in human_chromosomes() {
        let chr_str = version.versioned_chromosome(&chromosome);
        let Some(breakends) = map.get(&chr_str) else { continue };
        for b in breakends {
            let var = &variants[b.var as usize];
            if let Err(e) = output.write_breakend(&config, &genotype_ids, var, var.breakend(b.se as usize)) {
                failure = Some(e);
                break 'chromosomes;
            }
        }
    }
    if let Some(e) = failure {
        output.abort();
        return Err(e);
    }
    output.close()?;

    if config.write_breakend_tsv {
        write_breakend_tsv(&config, &variants, &map)?;
    }

    if crate::common::log::debug_enabled() {
        log_filter_type_counts(&variants);
    }

    info!("Esvee caller complete, mins({})", crate::common::fmt::java_fixed(start.elapsed().as_millis() as f64 / 60000.0, 3));
    Ok(())
}
