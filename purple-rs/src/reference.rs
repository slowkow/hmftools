//! ReferenceData: the reference genome version, the Ensembl data cache, the driver gene panel with the
//! DNDS likelihoods (DriverGeneResource), the somatic and germline hotspots, the cohort germline amp/del
//! frequencies and the target regions.

use anyhow::{bail, Result};
use std::collections::{BTreeMap, HashMap};

use hmf_common::config_builder::ConfigBuilder;
use hmf_common::driver_panel::DriverGene;
use hmf_common::ensembl::EnsemblDataCache;
use hmf_common::regions::RefGenomeVersion;
use hmf_common::variant_context::VcfReader;

use crate::config::PurpleConfig;
use crate::types::*;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DndsImpactLikelihood {
    pub drivers_per_sample: f64,
    pub passengers_per_mutation: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DndsGeneLikelihood {
    pub gene: String,
    pub missense: DndsImpactLikelihood,
    pub nonsense: DndsImpactLikelihood,
    pub splice: DndsImpactLikelihood,
    pub indel: DndsImpactLikelihood,
}

impl DndsGeneLikelihood {
    /// DndsDriverGeneLikelihood.select
    pub fn select(&self, impact: DriverImpact) -> DndsImpactLikelihood {
        match impact {
            DriverImpact::Nonsense => self.nonsense,
            DriverImpact::Missense => self.missense,
            DriverImpact::Inframe | DriverImpact::Frameshift => self.indel,
            DriverImpact::Splice => self.splice,
            DriverImpact::Unknown => panic!("IllegalArgumentException: no DNDS likelihood for impact UNKNOWN"),
        }
    }
}

/// DndsDriverGeneLikelihoodFile.fromInputStream
fn read_dnds(text: &str) -> Vec<DndsGeneLikelihood> {
    let p = |s: &str| hmf_common::java_text::parse_double(s).unwrap_or_else(|| panic!("NumberFormatException: {s}"));
    text.lines()
        .filter(|l| !l.starts_with("gene"))
        .map(|line| {
            let v: Vec<&str> = line.split('\t').collect();
            let imp = |o: usize| DndsImpactLikelihood { drivers_per_sample: p(v[o]), passengers_per_mutation: p(v[o + 1]) };
            DndsGeneLikelihood { gene: v[0].to_string(), missense: imp(1), nonsense: imp(3), splice: imp(5), indel: imp(7) }
        })
        .collect()
}

#[derive(Default)]
pub struct DriverGeneResource {
    pub driver_gene_list: Vec<DriverGene>,
    pub driver_gene_map: HashMap<String, usize>,
    pub tsg_likelihood: HashMap<String, DndsGeneLikelihood>,
    pub onco_likelihood: HashMap<String, DndsGeneLikelihood>,
}

impl DriverGeneResource {
    pub fn new(genes: Vec<DriverGene>) -> DriverGeneResource {
        let mut map = HashMap::new();
        for (i, g) in genes.iter().enumerate() {
            map.insert(g.gene.clone(), i);
        }
        let filter = |list: Vec<DndsGeneLikelihood>| -> HashMap<String, DndsGeneLikelihood> {
            let mut out = HashMap::new();
            for l in list.into_iter().filter(|l| map.contains_key(&l.gene)) {
                if out.contains_key(&l.gene) {
                    panic!("IllegalStateException: Duplicate key {}", l.gene);
                }
                out.insert(l.gene.clone(), l);
            }
            out
        };
        let tsg = filter(read_dnds(include_str!("../resources/dnds_driver_likelihood_tsg.tsv")));
        let onco = filter(read_dnds(include_str!("../resources/dnds_driver_likelihood_onco.tsv")));
        DriverGeneResource { driver_gene_list: genes, driver_gene_map: map, tsg_likelihood: tsg, onco_likelihood: onco }
    }
    pub fn gene(&self, name: &str) -> Option<&DriverGene> {
        self.driver_gene_map.get(name).map(|&i| &self.driver_gene_list[i])
    }
}

/// SimpleVariant of a hotspot file
#[derive(Clone, Debug, PartialEq)]
pub struct Hotspot {
    pub chromosome: String,
    pub position: i32,
    pub ref_bases: String,
    pub alt: String,
}

/// VariantHotspotFile.loadHotspotFile (the VCF form): by chromosome, in file order
pub fn load_hotspots(path: &str) -> Result<BTreeMap<Chromosome, Vec<Hotspot>>> {
    let mut out: BTreeMap<Chromosome, Vec<Hotspot>> = BTreeMap::new();
    if path.is_empty() {
        return Ok(out);
    }
    if !path.ends_with(".vcf.gz") {
        bail!("hotspot files other than .vcf.gz are not supported: {path}");
    }
    let mut reader = VcfReader::open(path)?;
    while let Some(r) = reader.next_record()? {
        if Chromosome::contains(&r.chrom) {
            out.entry(chr(&r.chrom)).or_default().push(Hotspot { chromosome: r.chrom.clone(), position: r.pos, ref_bases: r.ref_allele.clone(), alt: r.alt().to_string() });
        }
    }
    Ok(out)
}

/// AmpDelRegionFrequency
#[derive(Clone, Debug)]
pub struct AmpDelRegionFrequency {
    pub start: i32,
    pub end: i32,
    pub is_amp: bool,
    pub frequency: i32,
}

/// GermlineAmpDelFrequencyCache
#[derive(Default)]
pub struct GermlineAmpDelFrequencyCache {
    pub regions: HashMap<String, Vec<AmpDelRegionFrequency>>,
}

impl GermlineAmpDelFrequencyCache {
    pub fn load(path: Option<&str>) -> Result<GermlineAmpDelFrequencyCache> {
        let mut cache = GermlineAmpDelFrequencyCache::default();
        let Some(path) = path else { return Ok(cache) };
        let lines = hmf_common::java_text::read_lines(path)?;
        let Some(header) = lines.first() else { return Ok(cache) };
        let m = hmf_common::java_text::fields_index_map(header, ',');
        let idx = |k: &str| m.get(k).copied().ok_or_else(|| anyhow::anyhow!("{path}: missing column {k}"));
        let (ci, si, ei, fi) = (idx("Chromosome")?, idx("RegionStart")?, idx("RegionEnd")?, idx("Frequency")?);
        let ti = m.get("Type").copied();
        for line in &lines[1..] {
            let v: Vec<&str> = line.split(',').collect();
            let int = |s: &str| hmf_common::java_text::parse_int(s).ok_or_else(|| anyhow::anyhow!("For input string: \"{s}\""));
            let is_amp = match ti {
                Some(t) => match v[t] {
                    "AMP" => true,
                    "DEL" => false,
                    other => bail!("No enum constant EventType.{other}"),
                },
                None => false,
            };
            cache.regions.entry(v[ci].to_string()).or_default().push(AmpDelRegionFrequency { start: int(v[si])?, end: int(v[ei])?, is_amp, frequency: int(v[fi])? });
        }
        crate::info!("loaded {} germline amp-del frequencies from file: {}", lines.len() - 1, path);
        Ok(cache)
    }

    /// getRegionFrequency
    pub fn region_frequency(&self, chromosome: &str, start: i32, end: i32, buffer: i32, is_amp: bool) -> i32 {
        let Some(list) = self.regions.get(chromosome) else { return 0 };
        let mut total = 0;
        for r in list.iter().filter(|r| r.is_amp == is_amp) {
            if start - buffer > r.start {
                continue;
            }
            if r.start - buffer > start {
                break;
            }
            if (r.start - start).abs() <= buffer && (r.end - end).abs() <= buffer {
                total += r.frequency;
            }
        }
        total
    }
}

/// TargetRegionsData (the target regions and the TMB/TML settings)
#[derive(Default)]
pub struct TargetRegionsData {
    pub regions: HashMap<String, Vec<(i32, i32)>>,
    /// the regions with their BED name column (TaggedRegion), by versioned chromosome
    pub tagged: HashMap<String, Vec<(i32, i32, String)>>,
    pub total_bases: i32,
    pub chromosome_x_bases: i32,
    pub coding_bases: i32,
    pub coding_base_factor: i32,
    pub tmb_ratio: f64,
    pub tml_ratio: f64,
}

impl TargetRegionsData {
    pub fn has_target_regions(&self) -> bool {
        !self.regions.is_empty()
    }
    pub fn in_target_regions(&self, chromosome: &str, position: i32) -> bool {
        self.regions.get(chromosome).is_some_and(|l| l.iter().any(|(s, e)| position >= *s && position <= *e))
    }
    pub fn chromosome_x_region_percentage(&self) -> Option<f64> {
        if self.has_target_regions() { Some(self.chromosome_x_bases as f64 / self.total_bases as f64) } else { None }
    }

    /// TargetRegionsData.loadTargetRegionsBed
    pub fn load_bed(&mut self, path: &str, ensembl: &EnsemblDataCache) -> Result<()> {
        // TaggedRegion.loadRegionsFromBedFile: chromosome, start (0-based), end, optional tag; lines naming the
        // chromosome column skipped; then ChrBaseRegion.checkMergeOverlaps: sorted by chromosome rank and start,
        // regions within 1 base merged into the first
        let mut all: Vec<(Chromosome, i32, i32, String, i32)> = Vec::new();
        for line in hmf_common::java_text::read_lines(path)? {
            if line.contains("chromosome") || line.starts_with('#') || line.starts_with("track") || line.trim().is_empty() {
                continue;
            }
            let v: Vec<&str> = line.split('\t').collect();
            if v.len() < 3 || !Chromosome::contains(v[0]) {
                continue;
            }
            let start = hmf_common::java_text::parse_int(v[1]).unwrap_or(0) + 1;
            let end = hmf_common::java_text::parse_int(v[2]).unwrap_or(0);
            all.push((chr(v[0]), start, end, v.get(3).map(|t| t.to_string()).unwrap_or_default(), chromosome_rank(v[0])));
        }
        all.sort_by(|a, b| a.4.cmp(&b.4).then(a.1.cmp(&b.1)));
        let mut merged: Vec<(Chromosome, i32, i32, String, i32)> = Vec::new();
        for r in all {
            match merged.last_mut() {
                Some(last) if last.0 == r.0 && last.2 >= r.1 - 2 => {
                    if r.2 > last.2 {
                        last.2 = r.2;
                    }
                }
                _ => merged.push(r),
            }
        }
        let mut by_chr: BTreeMap<Chromosome, Vec<(i32, i32)>> = BTreeMap::new();
        let mut tagged: BTreeMap<Chromosome, Vec<(i32, i32, String)>> = BTreeMap::new();
        for (c, s, e, t, _) in merged {
            by_chr.entry(c).or_default().push((s, e));
            tagged.entry(c).or_default().push((s, e, t));
        }
        let version = ensembl.ref_genome_version();
        let excluded = crate::constants::TARGETED_TMB_GENE_EXCLUSIONS;
        let mut count = 0;
        for (c, regions) in &by_chr {
            let chr_str = version.versioned_chromosome(c.short_name());
            count += regions.len();
            let genes = ensembl.chr_gene_data_map().get(&chr_str).cloned().unwrap_or_default();
            let mut transcripts = Vec::new();
            for g in &genes {
                if excluded.contains(&g.gene_name.as_str()) {
                    continue;
                }
                if regions.iter().any(|(s, e)| positions_overlap(*s, *e, g.gene_start, g.gene_end)) {
                    if let Some(t) = ensembl.canonical_transcript_data(&g.gene_id) {
                        if !t.non_coding() {
                            transcripts.push(t.clone());
                        }
                    }
                }
            }
            for (s, e) in regions {
                self.total_bases += e - s + 1;
                if is_x(*c) {
                    self.chromosome_x_bases += e - s + 1;
                }
                let mut exonic: Vec<hmf_common::regions::BaseRegion> = Vec::new();
                for t in &transcripts {
                    let (cs, ce) = (t.coding_start.unwrap_or(0), t.coding_end.unwrap_or(0));
                    for exon in &t.exons {
                        if cs > exon.end {
                            continue;
                        }
                        if exon.start > ce {
                            break;
                        }
                        let max_cs = cs.max(exon.start);
                        let min_ce = ce.min(exon.end);
                        if positions_overlap(max_cs, min_ce, *s, *e) {
                            exonic.push(hmf_common::regions::BaseRegion::new(max_cs.max(*s), min_ce.min(*e)));
                        }
                    }
                }
                hmf_common::regions::merge_base_regions(&mut exonic, false);
                self.coding_bases += exonic.iter().map(|r| r.base_length()).sum::<i32>();
            }
            self.tagged.insert(chr_str.clone(), tagged.remove(c).unwrap_or_default());
            self.regions.insert(chr_str, regions.clone());
        }
        crate::info!("loaded {} target regions bases(total={} coding={}) from file({})", count, self.total_bases, self.coding_bases, path);
        Ok(())
    }
}

pub struct ReferenceData {
    pub version: RefGenomeVersion,
    pub v38: bool,
    pub ensembl: EnsemblDataCache,
    pub driver_genes: DriverGeneResource,
    pub other_reportable_transcripts: HashMap<String, Vec<String>>,
    pub somatic_hotspots: BTreeMap<Chromosome, Vec<Hotspot>>,
    pub germline_hotspots: BTreeMap<Chromosome, Vec<Hotspot>>,
    pub germline_amp_del_frequencies: GermlineAmpDelFrequencyCache,
    pub target_regions: TargetRegionsData,
    pub special_regions: crate::region::SpecialRegions,
}

impl ReferenceData {
    /// ReferenceData(configBuilder, config); the bool is isValid
    pub fn load(cb: &ConfigBuilder, config: &PurpleConfig) -> Result<(ReferenceData, bool)> {
        let mut valid = true;
        if !cb.has_value("ref_genome") {
            valid = false;
            crate::error!("ref_genome is a mandatory argument");
        }
        let version = if cb.has_value("ref_genome_version") {
            RefGenomeVersion::from_java(&cb.value("ref_genome_version").unwrap())?
        } else {
            // RefGenomeSource.deriveRefGenomeVersion: from the reference genome's chromosome names
            let fasta = cb.value("ref_genome").unwrap_or_default();
            let fai = std::fs::read_to_string(format!("{fasta}.fai")).unwrap_or_default();
            if fai.lines().any(|l| l.starts_with("chr1\t")) { RefGenomeVersion::from_java("38")? } else { RefGenomeVersion::from_java("37")? }
        };
        crate::info!("using ref genome: {}", version.name());
        let v38 = version.is_38();
        let somatic_hotspot_vcf = cb.value("somatic_hotspots").unwrap_or_default();
        let germline_hotspot_vcf = cb.value("germline_hotspots").unwrap_or_default();
        let mut driver_genes = DriverGeneResource::default();
        if config.run_drivers {
            if !config.germline_mode() {
                if somatic_hotspot_vcf.is_empty() {
                    valid = false;
                    crate::error!("somatic_hotspots is a mandatory argument when running drivers");
                }
                if !std::path::Path::new(&somatic_hotspot_vcf).exists() {
                    valid = false;
                    crate::error!("unable to open somatic_hotspots file {somatic_hotspot_vcf}");
                }
            }
            match hmf_common::driver_panel::read_driver_genes(&cb.value("driver_gene_panel").unwrap()) {
                Ok(genes) => driver_genes = DriverGeneResource::new(genes),
                Err(e) => {
                    valid = false;
                    crate::error!("unable to load driver genes: {e}");
                }
            }
            if cb.has_value("germline_vcf") {
                if germline_hotspot_vcf.is_empty() {
                    valid = false;
                    crate::error!("germline_hotspots is a mandatory argument when running drivers");
                }
                if !std::path::Path::new(&germline_hotspot_vcf).exists() {
                    valid = false;
                    crate::error!("unable to open germline_hotspots file {germline_hotspot_vcf}");
                }
            }
        }
        let mut other_reportable_transcripts = HashMap::new();
        for g in &driver_genes.driver_gene_list {
            if !g.additional_reported_transcripts.is_empty() {
                other_reportable_transcripts.insert(g.gene.clone(), g.additional_reported_transcripts.clone());
            }
        }
        let mut ensembl = EnsemblDataCache::new(&cb.value("ensembl_data_dir").unwrap(), version);
        ensembl.set_required_data(true, false, false, true);
        if !other_reportable_transcripts.is_empty() {
            let additional: Vec<String> = other_reportable_transcripts.values().flatten().cloned().collect();
            ensembl.load(true)?;
            ensembl.load_transcript_data(&[], &additional)?;
        } else {
            ensembl.load(false)?;
        }
        ensembl.create_gene_name_id_map();
        let invalid = hmf_common::driver_panel::find_invalid_driver_genes(&driver_genes.driver_gene_list, &ensembl);
        if !invalid.is_empty() {
            crate::error!("Invalid non-Ensembl gene names: {}", invalid.iter().map(|g| g.gene.as_str()).collect::<Vec<_>>().join(";"));
            std::process::exit(1);
        }
        let somatic_hotspots = match load_hotspots(&somatic_hotspot_vcf) {
            Ok(h) => h,
            Err(e) => {
                valid = false;
                crate::error!("failed to load hotspots: {e}");
                BTreeMap::new()
            }
        };
        let germline_hotspots = match load_hotspots(&germline_hotspot_vcf) {
            Ok(h) => h,
            Err(e) => {
                valid = false;
                crate::error!("failed to load hotspots: {e}");
                BTreeMap::new()
            }
        };
        let freq_file = if config.run_germline() { cb.value("germline_amp_del_freq_file") } else { None };
        let germline_amp_del_frequencies = GermlineAmpDelFrequencyCache::load(freq_file.as_deref())?;
        let mut target_regions = TargetRegionsData {
            coding_base_factor: cb.integer("target_coding_factor"),
            tmb_ratio: cb.decimal("target_tmb_ratio"),
            tml_ratio: cb.decimal("target_tml_ratio"),
            ..Default::default()
        };
        if let Some(bed) = cb.value("target_regions_bed") {
            target_regions.load_bed(&bed, &ensembl)?;
        }
        Ok((
            ReferenceData {
                version,
                v38,
                ensembl,
                driver_genes,
                other_reportable_transcripts,
                somatic_hotspots,
                germline_hotspots,
                germline_amp_del_frequencies,
                target_regions,
                special_regions: crate::region::SpecialRegions::new(v38),
            },
            valid,
        ))
    }
}
