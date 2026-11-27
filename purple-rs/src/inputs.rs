//! AmberData and CobaltData: the AMBER BAFs, QC and segments, the COBALT ratios, segments and the
//! CobaltChromosomes (gender and germline aberrations) derived from them.

use anyhow::{anyhow, bail, Context, Result};
use std::collections::BTreeMap;

use hmf_common::config_builder::check_add_dir_separator;
use hmf_common::java_text::{for_each_line, parse_double, parse_int};
use hmf_common::segmentation::pcf::{self, PcfPosition, PcfSource};

use crate::types::*;

pub const WINDOW_SIZE: i32 = 1000;

#[derive(Clone, Debug, PartialEq)]
pub struct AmberBaf {
    pub chromosome: String,
    pub position: i32,
    pub tumor_baf: f64,
    pub tumor_depth: i32,
    pub normal_baf: f64,
    pub normal_depth: i32,
}

impl AmberBaf {
    pub fn tumor_modified_baf(&self) -> f64 {
        0.5 + (self.tumor_baf - 0.5).abs()
    }
}

pub struct AmberData {
    /// by chromosome, in file order
    pub chromosome_bafs: BTreeMap<Chromosome, Vec<AmberBaf>>,
    /// the order in which Java's ArrayListMultimap first saw each chromosome (for the hash-order effects)
    pub chromosome_order: Vec<Chromosome>,
    pub tumor_segments: BTreeMap<Chromosome, Vec<PcfPosition>>,
    pub average_tumor_depth: i32,
    pub contamination: f64,
    pub patient_gender: Gender,
}

/// FileReaderUtils.getValue(lines, key, default, delim): the value after the first field equal to the key
fn qc_value(lines: &[String], key: &str) -> Option<String> {
    lines.iter().find_map(|l| {
        let mut parts = l.split('\t');
        if parts.next() == Some(key) { parts.next().map(|s| s.to_string()) } else { None }
    })
}

impl AmberData {
    pub fn load(sample: &str, amber_dir: &str, germline_only: bool, v38: bool, chr_x_region_percentage: Option<f64>) -> Result<AmberData> {
        let dir = check_add_dir_separator(amber_dir);
        let qc_file = format!("{dir}{sample}.amber.qc");
        if !std::path::Path::new(&qc_file).exists() {
            bail!("Unable to open Amber QC file: {qc_file}");
        }
        let lines = hmf_common::java_text::read_lines(&qc_file).with_context(|| qc_file.clone())?;
        let contamination = parse_double(&qc_value(&lines, "Contamination").unwrap_or_else(|| "0".into()))
            .ok_or_else(|| anyhow!("Unable to parse amber qc file with {} lines.", lines.len()))?;

        let baf_file = format!("{dir}{sample}.amber.baf.tsv.gz");
        if !std::path::Path::new(&baf_file).exists() {
            bail!("Unable to open Amber BAF file: {baf_file}");
        }
        let (chromosome_bafs, chromosome_order) = read_bafs(&baf_file, !germline_only)?;
        let tumor_segments = if !germline_only {
            let pcf_file = pcf::generate_baf_filename(&dir, sample);
            if !std::path::Path::new(&pcf_file).exists() {
                bail!("Unable to open Amber PCF file: {pcf_file}");
            }
            pcf::read_positions(WINDOW_SIZE, PcfSource::TumorBaf, &pcf_file)?
        } else {
            BTreeMap::new()
        };
        // IntStream.average of the positive depths (long sum), Math.round, default 100
        let mut sum: i64 = 0;
        let mut count: i64 = 0;
        for bafs in chromosome_bafs.values() {
            for b in bafs {
                if b.tumor_depth > 0 {
                    sum += b.tumor_depth as i64;
                    count += 1;
                }
            }
        }
        let average_tumor_depth = if count > 0 { jround(sum as f64 / count as f64) as i32 } else { 100 };
        let patient_gender = determine_gender(v38, &chromosome_bafs, chr_x_region_percentage);
        Ok(AmberData { chromosome_bafs, chromosome_order, tumor_segments, average_tumor_depth, contamination, patient_gender })
    }

    pub fn min_somatic_total_read_count(&self) -> i32 {
        (0.8 * self.average_tumor_depth as f64).floor() as i32
    }
    pub fn max_somatic_total_read_count(&self) -> i32 {
        (1.2 * self.average_tumor_depth as f64).ceil() as i32
    }

    #[cfg(test)]
    pub fn for_test(average_tumor_depth: i32, gender: Gender) -> AmberData {
        AmberData {
            chromosome_bafs: BTreeMap::new(),
            chromosome_order: Vec::new(),
            tumor_segments: BTreeMap::new(),
            average_tumor_depth,
            contamination: 0.0,
            patient_gender: gender,
        }
    }
}

fn column(map: &std::collections::HashMap<String, usize>, name: &str, file: &str) -> Result<usize> {
    map.get(name).copied().ok_or_else(|| anyhow!("{file}: missing column {name}"))
}

/// AmberBAFFile.read
fn read_bafs(file: &str, has_tumor: bool) -> Result<(BTreeMap<Chromosome, Vec<AmberBaf>>, Vec<Chromosome>)> {
    let mut map: BTreeMap<Chromosome, Vec<AmberBaf>> = BTreeMap::new();
    let mut order: Vec<Chromosome> = Vec::new();
    let mut idx: Option<[usize; 6]> = None;
    for_each_line(file, |line| {
        let Some(ix) = idx else {
            let m = hmf_common::java_text::fields_index_map(line, '\t');
            idx = Some([
                column(&m, "chromosome", file)?,
                column(&m, "position", file)?,
                column(&m, "tumorBAF", file)?,
                column(&m, "tumorDepth", file)?,
                column(&m, "normalBAF", file)?,
                column(&m, "normalDepth", file)?,
            ]);
            return Ok(());
        };
        let v: Vec<&str> = line.split('\t').collect();
        let num = |i: usize| -> Result<f64> { parse_double(v[i]).ok_or_else(|| anyhow!("For input string: \"{}\"", v[i])) };
        let int = |i: usize| -> Result<i32> { parse_int(v[i]).ok_or_else(|| anyhow!("For input string: \"{}\"", v[i])) };
        let chromosome = v[ix[0]];
        let c = Chromosome::parse(chromosome)?;
        let baf = AmberBaf {
            chromosome: chromosome.to_string(),
            position: int(ix[1])?,
            tumor_baf: if has_tumor { num(ix[2])? } else { 0.5 },
            tumor_depth: int(ix[3])?,
            normal_baf: num(ix[4])?,
            normal_depth: int(ix[5])?,
        };
        if !map.contains_key(&c) {
            order.push(c);
        }
        map.entry(c).or_default().push(baf);
        Ok(())
    })?;
    Ok((map, order))
}

/// AmberGender.determineGender
fn determine_gender(v38: bool, bafs: &BTreeMap<Chromosome, Vec<AmberBaf>>, chr_x_region_percentage: Option<f64>) -> Gender {
    let (start, end) = if v38 { (2781479, 156030895) } else { (2699520, 155260560) };
    let mut total: i64 = 0;
    let mut inclusion: i64 = 0;
    for (c, list) in bafs {
        total += list.len() as i64;
        if is_x(*c) {
            inclusion = list.iter().filter(|b| b.position >= start && b.position <= end).count() as i64;
        }
    }
    if inclusion < 3 {
        return Gender::Male;
    }
    let perc = inclusion as f64 / total as i32 as f64;
    let threshold = match chr_x_region_percentage {
        Some(p) => 0.05 * p,
        None => 0.01,
    };
    if perc > threshold { Gender::Female } else { Gender::Male }
}

// ---------------------------------------------------------------------------------------------
// COBALT

#[derive(Clone, Debug, PartialEq)]
pub struct CobaltRatio {
    pub chromosome: String,
    pub position: i32,
    pub reference_read_depth: f64,
    pub reference_gc_ratio: f64,
    pub reference_gc_content: f64,
    pub reference_gc_diploid_ratio: f64,
    pub tumor_read_depth: f64,
    pub tumor_gc_ratio: f64,
    pub tumor_gc_content: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CobaltChromosome {
    pub chromosome: Chromosome,
    pub typical_ratio: f64,
    pub actual_ratio: f64,
    pub mosiac: bool,
}

impl CobaltChromosome {
    pub fn is_normal(&self) -> bool {
        equal(self.typical_ratio, self.actual_ratio)
    }
    pub fn is_diploid(&self) -> bool {
        equal(self.actual_ratio, 1.0)
    }
    pub fn is_autosome(&self) -> bool {
        is_autosome(self.chromosome)
    }
}

#[derive(Clone, Debug)]
pub struct MedianRatio {
    pub chromosome: String,
    pub median_ratio: f64,
    pub count: i32,
}

pub struct CobaltChromosomes {
    pub gender: Gender,
    /// the aberrations found, in the order of GermlineAberration (Java: a HashSet of the enum)
    pub aberrations: Vec<GermlineAberration>,
    /// by contig name as in the ratio file
    pub chromosomes: std::collections::HashMap<String, CobaltChromosome>,
    /// the same entries by chromosome ordinal, with their contig name
    by_ordinal: Vec<Option<(String, CobaltChromosome)>>,
}

impl CobaltChromosomes {
    /// CobaltChromosomes(medianRatios, calcAberrations)
    pub fn new(median_ratios: &[MedianRatio], calc_aberrations: bool) -> CobaltChromosomes {
        let mut min_autosome_ratio = 0.0;
        let mut y_median = 0.0;
        let mut x_median = 0.0;
        for m in median_ratios {
            let c = chr(&m.chromosome);
            if is_autosome(c) {
                min_autosome_ratio = if min_autosome_ratio == 0.0 { m.median_ratio } else { jmin(min_autosome_ratio, m.median_ratio) };
            } else if is_y(c) && m.count >= 1000 {
                y_median = m.median_ratio;
            } else if is_x(c) {
                x_median = m.median_ratio;
            }
        }
        let is_female = greater_or_equal(x_median, 0.65) && less_than(y_median, 0.05);
        let gender = if is_female { Gender::Female } else { Gender::Male };
        let mut aberrations: Vec<GermlineAberration> = Vec::new();
        let mut chromosomes = std::collections::HashMap::new();
        for m in median_ratios {
            let c = chr(&m.chromosome);
            if is_y(c) && m.count < 1000 {
                continue;
            }
            let aberration = if calc_aberrations { aberration(is_female, c, m.median_ratio, min_autosome_ratio) } else { GermlineAberration::None };
            let typical = if is_female {
                if is_y(c) { 0.0 } else { 1.0 }
            } else if is_autosome(c) {
                1.0
            } else {
                0.5
            };
            let actual = match aberration {
                GermlineAberration::TrisomyX
                | GermlineAberration::Trisomy13
                | GermlineAberration::Trisomy15
                | GermlineAberration::Trisomy18
                | GermlineAberration::Trisomy21 => 1.5,
                GermlineAberration::Tertrasomy9 => 2.0,
                GermlineAberration::MosaicX => m.median_ratio,
                GermlineAberration::Xyy | GermlineAberration::Klinefelter => 1.0,
                GermlineAberration::None => typical,
            };
            if aberration != GermlineAberration::None && !aberrations.contains(&aberration) {
                aberrations.push(aberration);
            }
            if positive(typical) {
                chromosomes.insert(
                    m.chromosome.clone(),
                    CobaltChromosome { chromosome: c, typical_ratio: typical, actual_ratio: actual, mosiac: aberration == GermlineAberration::MosaicX },
                );
            }
        }
        if aberrations.is_empty() {
            aberrations.push(GermlineAberration::None);
        }
        aberrations.sort();
        let mut by_ordinal = vec![None; 24];
        for (k, v) in &chromosomes {
            if let Some(c) = Chromosome::from_name(k) {
                by_ordinal[c.ordinal()] = Some((k.clone(), v.clone()));
            }
        }
        CobaltChromosomes { gender, aberrations, chromosomes, by_ordinal }
    }

    pub fn has_chromosome(&self, contig: &str) -> bool {
        self.get(contig).is_some()
    }
    /// by the chromosome's ordinal (the map is keyed by the COBALT contig names, which the regions share)
    pub fn get(&self, contig: &str) -> Option<&CobaltChromosome> {
        let c = Chromosome::from_name(contig)?;
        let found = self.by_ordinal[c.ordinal()].as_ref()?;
        if found.0 == contig { Some(&found.1) } else { self.chromosomes.get(contig) }
    }
}

fn aberration(is_female: bool, c: Chromosome, median: f64, min_autosome_ratio: f64) -> GermlineAberration {
    let name = c.short_name();
    if greater_or_equal(median, 1.35) && matches!(name, "13" | "15" | "18" | "21" | "X") {
        return match name {
            "13" => GermlineAberration::Trisomy13,
            "15" => GermlineAberration::Trisomy15,
            "18" => GermlineAberration::Trisomy18,
            "21" => GermlineAberration::Trisomy21,
            _ => GermlineAberration::TrisomyX,
        };
    }
    if !is_female && is_y(c) && greater_or_equal(median, 0.75) {
        return GermlineAberration::Xyy;
    }
    if is_female && less_than(median, jmin(min_autosome_ratio, 0.8)) && is_x(c) {
        return GermlineAberration::MosaicX;
    }
    if !is_female && is_x(c) && greater_or_equal(median, 0.65) {
        return GermlineAberration::Klinefelter;
    }
    GermlineAberration::None
}

pub struct CobaltData {
    pub cobalt_chromosomes: CobaltChromosomes,
    pub ratios: BTreeMap<Chromosome, Vec<CobaltRatio>>,
    pub tumor_segments: BTreeMap<Chromosome, Vec<PcfPosition>>,
    pub reference_segments: BTreeMap<Chromosome, Vec<PcfPosition>>,
}

impl CobaltData {
    pub fn gender(&self) -> Gender {
        self.cobalt_chromosomes.gender
    }

    pub fn load(
        reference_id: Option<&str>,
        tumor_id: Option<&str>,
        cobalt_dir: &str,
        amber_gender: Gender,
        tumor_only: bool,
        germline_only: bool,
    ) -> Result<CobaltData> {
        let dir = check_add_dir_separator(cobalt_dir);
        let tumor_segments = match tumor_id {
            Some(t) => pcf::read_positions(WINDOW_SIZE, PcfSource::TumorRatio, &pcf::generate_ratio_filename(&dir, t))?,
            None => BTreeMap::new(),
        };
        let reference_segments = match reference_id {
            Some(r) => pcf::read_positions(WINDOW_SIZE, PcfSource::ReferenceRatio, &pcf::generate_ratio_filename(&dir, r))?,
            None => BTreeMap::new(),
        };
        let tumor = tumor_id.unwrap_or("null");
        let mut file = format!("{dir}{tumor}.cobalt.ratio.tsv.gz");
        if !std::path::Path::new(&file).exists() {
            let unzipped = file.replace(".gz", "");
            if std::path::Path::new(&unzipped).exists() {
                file = unzipped;
            }
        }
        if !std::path::Path::new(&file).exists() {
            bail!("unable to open Cobalt ratio file: {file}");
        }
        let ratios = read_ratios(&file, if tumor_only { Some(amber_gender) } else { None }, !germline_only)?;
        let median_ratios = median_ratios(&ratios);
        let cobalt_chromosomes = CobaltChromosomes::new(&median_ratios, !tumor_only);
        Ok(CobaltData { cobalt_chromosomes, ratios, tumor_segments, reference_segments })
    }
}

/// CobaltRatioFile.genderAdjustedDiploidRatio
fn gender_adjusted(gender: Option<Gender>, contig: &str, ratio: f64) -> f64 {
    let Some(g) = gender else { return ratio };
    if !Chromosome::contains(contig) {
        return ratio;
    }
    let c = chr(contig);
    if is_x(c) {
        return if g == Gender::Female { 1.0 } else { 0.5 };
    }
    if is_y(c) {
        return if g == Gender::Female { 0.0 } else { 0.5 };
    }
    ratio
}

/// CobaltRatioFile.read(filename, gender, hasTumor)
fn read_ratios(file: &str, gender: Option<Gender>, has_tumor: bool) -> Result<BTreeMap<Chromosome, Vec<CobaltRatio>>> {
    let mut map: BTreeMap<Chromosome, Vec<CobaltRatio>> = BTreeMap::new();
    let mut idx: Option<[usize; 9]> = None;
    let mut current: Option<(String, Chromosome)> = None;
    for_each_line(file, |line| {
        let Some(ix) = idx else {
            let m = hmf_common::java_text::fields_index_map(line, '\t');
            idx = Some([
                column(&m, "chromosome", file)?,
                column(&m, "position", file)?,
                column(&m, "referenceReadDepth", file)?,
                column(&m, "tumorReadDepth", file)?,
                column(&m, "referenceGCRatio", file)?,
                column(&m, "tumorGCRatio", file)?,
                column(&m, "referenceGCDiploidRatio", file)?,
                column(&m, "referenceGCContent", file)?,
                column(&m, "tumorGCContent", file)?,
            ]);
            return Ok(());
        };
        let v: Vec<&str> = line.split('\t').collect();
        let num = |i: usize| -> Result<f64> { parse_double(v[i]).ok_or_else(|| anyhow!("For input string: \"{}\"", v[i])) };
        let chromosome = v[ix[0]];
        if current.as_ref().is_none_or(|(name, _)| name != chromosome) {
            let c = Chromosome::parse(chromosome)?;
            // a chromosome seen again starts a new list (HashMap.put)
            map.insert(c, Vec::new());
            current = Some((chromosome.to_string(), c));
        }
        let (name, c) = current.clone().unwrap();
        let ref_read_depth = num(ix[2])?;
        let mut ref_gc_ratio = num(ix[4])?;
        let mut ref_gc_diploid_ratio = num(ix[6])?;
        if ref_read_depth == -1.0 {
            ref_gc_ratio = 1.0;
            ref_gc_diploid_ratio = 1.0;
        }
        let ref_gc_ratio = gender_adjusted(gender, &name, ref_gc_ratio);
        let ref_gc_diploid_ratio = gender_adjusted(gender, &name, ref_gc_diploid_ratio);
        let tumor_gc_ratio = if has_tumor { num(ix[5])? } else { ref_gc_diploid_ratio };
        let ratio = CobaltRatio {
            chromosome: name,
            position: parse_int(v[ix[1]]).ok_or_else(|| anyhow!("For input string: \"{}\"", v[ix[1]]))?,
            reference_read_depth: ref_read_depth,
            reference_gc_ratio: ref_gc_ratio,
            reference_gc_content: num(ix[7])?,
            reference_gc_diploid_ratio: ref_gc_diploid_ratio,
            tumor_read_depth: num(ix[3])?,
            tumor_gc_ratio,
            tumor_gc_content: num(ix[8])?,
        };
        map.get_mut(&c).unwrap().push(ratio);
        Ok(())
    })?;
    Ok(map)
}

/// MedianRatioFactory.create(Map): the median of the positive reference GC ratios of each chromosome
/// (Java iterates its HashMap keyed by HumanChromosome; CobaltChromosomes does not depend on the order
/// except through which aberrations it collects, which are kept sorted here)
pub fn median_ratios(ratios: &BTreeMap<Chromosome, Vec<CobaltRatio>>) -> Vec<MedianRatio> {
    ratios
        .iter()
        .map(|(c, list)| {
            let name = list.first().map(|r| r.chromosome.clone()).unwrap_or_else(|| c.short_name().to_string());
            let positive_ratios: Vec<f64> = list.iter().map(|r| r.reference_gc_ratio).filter(|r| positive(*r)).collect();
            let count = positive_ratios.len() as i32;
            MedianRatio { chromosome: name, median_ratio: if count > 0 { median(&positive_ratios) } else { 0.0 }, count }
        })
        .collect()
}
