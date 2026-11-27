//! The haplotype panel: events (VariantHaplotypeEvent, HaplotypeEventFactory), haplotypes (DefaultHaplotype,
//! NonDefaultHaplotype), GeneHaplotypePanel and HaplotypePanel, and its loader (PanelLoader).

use std::collections::{BTreeMap, BTreeSet, HashSet};

use anyhow::{anyhow, bail, Context, Result};

/// HumanChromosome: "1".."22", "X", "Y" by index 0..24
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Chromosome(pub u8);

pub const CHROMOSOME_NAMES: [&str; 24] = hmf_common::regions::HUMAN_CHROMOSOMES;

impl Chromosome {
    /// HumanChromosome.fromString: a "chr" prefix of any case is dropped, the rest must be a HumanChromosome
    /// name exactly (valueOf("_" + name))
    pub fn from_string(contig: &str) -> Result<Chromosome> {
        let bare = if contig.len() >= 3 && contig.is_char_boundary(3) && contig[..3].eq_ignore_ascii_case("chr") { &contig[3..] } else { contig };
        CHROMOSOME_NAMES
            .iter()
            .position(|n| *n == bare)
            .map(|i| Chromosome(i as u8))
            .ok_or_else(|| anyhow!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._{bare}"))
    }
    pub fn name(self) -> &'static str {
        CHROMOSOME_NAMES[self.0 as usize]
    }
}

pub const EVENT_TYPE_STRING: &str = "VAR";
pub const ID_FIELD_COUNT: usize = 5;

/// VariantHaplotypeEvent, the only kind of HaplotypeEvent
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VariantHaplotypeEvent {
    pub chromosome: Chromosome,
    pub position: i32,
    pub ref_: String,
    pub alt: String,
}

impl VariantHaplotypeEvent {
    pub fn new(chromosome: Chromosome, position: i32, ref_: &str, alt: &str) -> VariantHaplotypeEvent {
        VariantHaplotypeEvent { chromosome, position, ref_: ref_.to_string(), alt: alt.to_string() }
    }

    /// "VAR_chr<chromosome>_<position>_<ref>_<alt>"
    pub fn id(&self) -> String {
        format!("{EVENT_TYPE_STRING}_chr{}_{}_{}_{}", self.chromosome.name(), self.position, self.ref_, self.alt)
    }

    /// HaplotypeEventFactory.fromId / VariantHaplotypeEvent.fromId: five fields separated by '_', the id
    /// written back the same
    pub fn from_id(event_id: &str) -> Result<VariantHaplotypeEvent> {
        let split = hmf_common::java_text::split(event_id, '_');
        if split[0] != EVENT_TYPE_STRING {
            bail!("java.lang.IllegalArgumentException: Cannot construct HaplotypeEvent from id: {event_id}");
        }
        if split.len() != ID_FIELD_COUNT {
            bail!(
                "java.lang.IllegalArgumentException: ID '{event_id}' of VariantHaplotypeEvent has incorrect field count: {} instead of {ID_FIELD_COUNT}",
                split.len()
            );
        }
        let chromosome = Chromosome::from_string(split[1])?;
        let position = hmf_common::java_text::parse_int(split[2])
            .ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{}\"", split[2]))?;
        let event = VariantHaplotypeEvent::new(chromosome, position, split[3], split[4]);
        if event.id() != event_id {
            bail!("java.lang.IllegalArgumentException: VariantHaplotypeEvent derived from event ID '{event_id}' has different ID '{}'", event.id());
        }
        Ok(event)
    }

    /// getCoveredPositions: position up to the last reference base
    pub fn covered_positions(&self) -> std::ops::Range<i32> {
        self.position..self.position + self.ref_.chars().count() as i32
    }

    /// isRelevantFor: on the same chromosome with overlapping reference bases
    pub fn is_relevant_for(&self, other: &VariantHaplotypeEvent) -> bool {
        other.chromosome == self.chromosome
            && self.position < other.position + other.ref_.chars().count() as i32
            && other.position < self.position + self.ref_.chars().count() as i32
    }
}

#[derive(Clone, Debug)]
pub struct DefaultHaplotype {
    pub name: String,
    pub is_wild_type: bool,
    pub events_to_ignore: Vec<VariantHaplotypeEvent>,
}

#[derive(Clone, Debug)]
pub struct NonDefaultHaplotype {
    pub name: String,
    pub is_wild_type: bool,
    pub events: Vec<VariantHaplotypeEvent>,
}

impl NonDefaultHaplotype {
    pub fn new(name: &str, is_wild_type: bool, events: Vec<VariantHaplotypeEvent>) -> Result<NonDefaultHaplotype> {
        if events.is_empty() {
            bail!("java.lang.RuntimeException: Non-wild type haplotype '{name}' has no associated events");
        }
        Ok(NonDefaultHaplotype { name: name.to_string(), is_wild_type, events })
    }
    pub fn is_relevant_for(&self, event: &VariantHaplotypeEvent) -> bool {
        self.events.iter().any(|e| event.is_relevant_for(e))
    }
    /// getMatchingCount: how many of the haplotype's events have this id
    pub fn matching_count(&self, event_id: &str) -> i32 {
        self.events.iter().filter(|e| e.id() == event_id).count() as i32
    }
}

#[derive(Clone, Debug)]
pub struct GeneHaplotypePanel {
    pub default_haplotype: DefaultHaplotype,
    pub non_default_haplotypes: Vec<NonDefaultHaplotype>,
    pub wild_type_haplotype_name: String,
}

impl GeneHaplotypePanel {
    /// isRelevantFor: not one of the default haplotype's events to ignore, and overlapping an event of a
    /// non-default haplotype
    pub fn is_relevant_for(&self, event: &VariantHaplotypeEvent) -> bool {
        let id = event.id();
        let ignored = self.default_haplotype.events_to_ignore.iter().any(|e| e.id() == id);
        !ignored && self.non_default_haplotypes.iter().any(|h| h.is_relevant_for(event))
    }
    pub fn haplotype_count(&self) -> usize {
        self.non_default_haplotypes.len() + 1
    }
}

/// HaplotypePanel: gene to its panel (Java's HashMap; nothing written depends on its order)
#[derive(Clone, Debug, Default)]
pub struct HaplotypePanel {
    pub genes: BTreeMap<String, GeneHaplotypePanel>,
}

impl HaplotypePanel {
    /// getRelevantVariantPositions: the reference positions of the non-default haplotypes' events
    pub fn relevant_variant_positions(&self) -> BTreeMap<Chromosome, HashSet<i32>> {
        let mut out: BTreeMap<Chromosome, HashSet<i32>> = BTreeMap::new();
        for panel in self.genes.values() {
            for h in &panel.non_default_haplotypes {
                for e in &h.events {
                    out.entry(e.chromosome).or_default().extend(e.covered_positions());
                }
            }
        }
        out
    }
    /// isRelevantFor(eventId, gene): the id parsed again as HaplotypeEventFactory.fromId does
    pub fn is_relevant_for_id(&self, event_id: &str, gene: &str) -> Result<bool> {
        let event = VariantHaplotypeEvent::from_id(event_id)?;
        Ok(self.genes[gene].is_relevant_for(&event))
    }
    pub fn haplotype_count(&self) -> usize {
        self.genes.values().map(|g| g.haplotype_count()).sum()
    }
}

pub const HAPLOTYPE_EVENT_DELIMITER: char = ';';

fn haplotype_events(text: &str) -> Result<Vec<VariantHaplotypeEvent>> {
    if text.is_empty() {
        return Ok(Vec::new());
    }
    hmf_common::java_text::split(text, HAPLOTYPE_EVENT_DELIMITER).into_iter().map(VariantHaplotypeEvent::from_id).collect()
}

fn column(fields: &std::collections::HashMap<String, usize>, name: &str) -> Result<usize> {
    fields.get(name).copied().ok_or_else(|| anyhow!("java.lang.NullPointerException: no column {name}"))
}

/// PanelLoader.loadHaplotypePanel: one line per haplotype (gene, haplotype, default, wildType and its events
/// separated by ';'); each gene has exactly one default haplotype and one wild type
pub fn load_haplotype_panel(filename: &str) -> Result<HaplotypePanel> {
    let r = (|| -> Result<HaplotypePanel> {
        let lines = hmf_common::java_text::read_lines(filename)?;
        let header = lines.first().ok_or_else(|| anyhow!("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0"))?;
        let fields = hmf_common::java_text::fields_index_map(header, '\t');
        let gene_i = column(&fields, "gene")?;
        let haplotype_i = column(&fields, "haplotype")?;
        let default_i = column(&fields, "default")?;
        let wild_type_i = column(&fields, "wildType")?;
        let events_i = column(&fields, "events")?;

        let mut defaults: BTreeMap<String, Vec<DefaultHaplotype>> = BTreeMap::new();
        let mut non_defaults: BTreeMap<String, Vec<NonDefaultHaplotype>> = BTreeMap::new();
        for line in &lines[1..] {
            let values = hmf_common::java_text::split_all(line, '\t');
            let get = |i: usize| -> Result<&str> {
                values.get(i).copied().ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {}", values.len()))
            };
            let events_text = get(events_i)?;
            let gene = get(gene_i)?.to_string();
            let is_default = hmf_common::java_text::parse_bool(get(default_i)?);
            let is_wild_type = hmf_common::java_text::parse_bool(get(wild_type_i)?);
            let events = haplotype_events(events_text)?;
            let name = get(haplotype_i)?;
            if is_default {
                defaults.entry(gene).or_default().push(DefaultHaplotype { name: name.to_string(), is_wild_type, events_to_ignore: events });
            } else {
                non_defaults.entry(gene).or_default().push(NonDefaultHaplotype::new(name, is_wild_type, events)?);
            }
        }

        let mut panel = HaplotypePanel::default();
        let genes: BTreeSet<String> = defaults.keys().chain(non_defaults.keys()).cloned().collect();
        for gene in genes {
            let gene_defaults = defaults.remove(&gene).unwrap_or_default();
            if gene_defaults.len() != 1 {
                bail!("java.lang.RuntimeException: Need to have exactly 1 default haplotype for gene: {gene}");
            }
            let default_haplotype = gene_defaults.into_iter().next().unwrap();
            let gene_non_defaults = non_defaults.remove(&gene).unwrap_or_default();
            if gene_non_defaults.iter().any(|h| h.name == default_haplotype.name) {
                bail!("java.lang.RuntimeException: Cannot have non-default haplotype with same name as default haplotype for gene: {gene}");
            }
            let wild_type_haplotype_name = wild_type_haplotype_name(&default_haplotype, &gene_non_defaults, &gene)?;
            panel.genes.insert(gene, GeneHaplotypePanel { default_haplotype, non_default_haplotypes: gene_non_defaults, wild_type_haplotype_name });
        }
        Ok(panel)
    })();
    let panel = r.with_context(|| format!("java.lang.RuntimeException: failed to load haplotypes TSV: {filename}"))?;
    hmf_common::info!("loaded {} haplotypes from file ({})", panel.haplotype_count(), filename);
    Ok(panel)
}

fn wild_type_haplotype_name(default: &DefaultHaplotype, non_defaults: &[NonDefaultHaplotype], gene: &str) -> Result<String> {
    let mut names: Vec<String> = non_defaults.iter().filter(|h| h.is_wild_type).map(|h| h.name.clone()).collect();
    if default.is_wild_type {
        names.push(default.name.clone());
    }
    if names.len() != 1 {
        bail!("java.lang.RuntimeException: Need to have exactly 1 wild type haplotype for gene: {gene}");
    }
    let name = names.remove(0);
    let clash = !default.is_wild_type && non_defaults.iter().any(|h| h.name == name && !h.is_wild_type);
    if clash {
        bail!("java.lang.RuntimeException: Wild type status inconsistent for haplotype {name} for gene {gene}");
    }
    Ok(name)
}
