//! IndelPrep, IndelVariant, IndelDetails and IndelContext: INDELs counted by type (deletion or insertion),
//! context (in a repeat, with microhomology, or neither) and length (or homology length), 30 bins.

use anyhow::{anyhow, bail, Result};

use hmf_common::debug;
use hmf_common::variant_context::{VariantType, VcfRecord};

use crate::compat;
use crate::prep::MutContextCount;
use crate::prep::LazyRefGenome;
use crate::refgenome::RefGenomeSource;
use crate::vcf_file::{type_name, SmallVariant};

pub const INDEL_DETAILS_FILE_SUFFIX: &str = ".chord.indel.details.tsv";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MutationType {
    Deletion,
    Insertion,
}

impl MutationType {
    pub const VALUES: [MutationType; 2] = [MutationType::Deletion, MutationType::Insertion];
    pub fn abbreviation(self) -> &'static str {
        match self {
            MutationType::Deletion => "del",
            MutationType::Insertion => "ins",
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            MutationType::Deletion => "DELETION",
            MutationType::Insertion => "INSERTION",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ContextType {
    Repeat,
    Microhomology,
    None,
}

impl ContextType {
    pub const VALUES: [ContextType; 3] = [ContextType::Repeat, ContextType::Microhomology, ContextType::None];
    pub fn suffix(self) -> &'static str {
        match self {
            ContextType::Repeat => ".rep.len.",
            ContextType::Microhomology => ".mh.bimh.",
            ContextType::None => ".none.len.",
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            ContextType::Repeat => "REPEAT",
            ContextType::Microhomology => "MICROHOMOLOGY",
            ContextType::None => "NONE",
        }
    }
}

pub const INDEL_LENGTH_CAP: i32 = 5;
pub const HOM_BASES_COUNT_CAP: i32 = 5;

// IndelDetails: thresholds for the context of an indel
const REPETITIVE_REGION_MIN_REPEAT_UNITS: i32 = 2;
const SHORT_REPEAT_UNIT_MAX_LENGTH: i32 = 50;
const MIN_HOMOLOGOUS_BASES: i32 = 2;
const LONG_INDEL_MIN_HOMOLOGOUS_BASES: i32 = 1;
const LONG_INDEL_MIN_LENGTH: i32 = 4;
const LEFT_FLANK_INDEL_LENGTH_UNITS: i32 = 1;
const RIGHT_FLANK_INDEL_LENGTH_UNITS: i32 = REPETITIVE_REGION_MIN_REPEAT_UNITS + 1;

/// IndelVariant
#[derive(Clone, Debug)]
pub struct IndelVariant {
    pub variant: SmallVariant,
    pub indel_length: i32,
    pub is_deletion: bool,
}

fn null_bases() -> anyhow::Error {
    anyhow!("java.lang.NullPointerException: Cannot invoke \"java.lang.String.<init>(byte[])\" because the bases are null")
}

impl IndelVariant {
    pub fn new(variant: SmallVariant) -> IndelVariant {
        let signed = variant.alt_bases.len() as i32 - variant.ref_bases.len() as i32;
        IndelVariant { variant, indel_length: signed.abs(), is_deletion: signed < 0 }
    }

    /// getFlankingRefBasesRight: the bases after the indel, `units` indel lengths long, cut at the contig end
    pub fn flanking_ref_bases_right(&self, ref_genome: &RefGenomeSource, units: i32) -> Result<String> {
        let pos = self.variant.position;
        let len = self.indel_length;
        let (start, end) = if self.is_deletion { (pos + 1 + len, pos + len + len * units) } else { (pos + 1, pos + len * units) };
        let end = end.min(ref_genome.chromosome_length(&self.variant.chromosome)?);
        match ref_genome.get_bases(&self.variant.chromosome, start, end)? {
            Some(b) => Ok(String::from_utf8_lossy(&b).into_owned()),
            // an indel at the contig end: no bases after it (CH5)
            None if !compat::java("CH5") && start > end => Ok(String::new()),
            None => Err(null_bases()),
        }
    }

    /// getFlankingRefBasesLeft: the bases up to and including the anchor base, `units` indel lengths long
    pub fn flanking_ref_bases_left(&self, ref_genome: &RefGenomeSource, units: i32) -> Result<String> {
        let mut start = self.variant.position - self.indel_length * units + 1;
        let end = self.variant.position;
        if !compat::java("CH5") {
            // an indel near the contig start: the bases from the first (CH5)
            start = start.max(1);
        }
        match ref_genome.get_bases(&self.variant.chromosome, start, end)? {
            Some(b) => Ok(String::from_utf8_lossy(&b).into_owned()),
            None => Err(null_bases()),
        }
    }

    /// getIndelSequence: the allele less its first base
    pub fn indel_sequence(&self) -> String {
        let s = if self.is_deletion { &self.variant.ref_bases } else { &self.variant.alt_bases };
        s.get(1..).unwrap_or("").to_string()
    }
}

/// IndelVariant.countHomologousBases: the common prefix of the two
pub fn count_homologous_bases(indel_sequence: &str, flank_sequence: &str) -> i32 {
    indel_sequence.bytes().zip(flank_sequence.bytes()).take_while(|(a, b)| a == b).count() as i32
}

/// IndelVariant.countRepeatUnits: how many copies of the indel sequence start the indel sequence plus
/// the right flank (Java replaces each match of the sequence, as a regex, by '#' and counts the leading
/// '#'s; indel sequences are bases, which have no regex meaning)
pub fn count_repeat_units(indel_sequence: &str, right_flank: &str) -> i32 {
    let sequence = format!("{indel_sequence}{right_flank}");
    if indel_sequence.is_empty() {
        // Java: the empty pattern matches before every character, so the sequence starts with one '#'
        return 1;
    }
    let mut count = 0;
    let mut rest = sequence.as_str();
    while let Some(r) = rest.strip_prefix(indel_sequence) {
        count += 1;
        rest = r;
    }
    count
}

#[derive(Clone, Debug)]
pub struct IndelDetails {
    pub chromosome: String,
    pub position: i32,
    pub ref_bases: String,
    pub alt_bases: String,
    pub indel_sequence: String,
    pub indel_length: i32,
    pub right_flank_bases: String,
    pub left_flank_bases: String,
    pub right_hom_bases_count: i32,
    pub left_hom_bases_count: i32,
    pub max_hom_bases_count: i32,
    pub repeat_units_count: i32,
    pub mutation_type: MutationType,
    pub context_type: ContextType,
}

fn reverse(s: &str) -> String {
    s.chars().rev().collect()
}

/// IndelDetails.determineIndelContext
fn determine_indel_context(indel_length: i32, repeat_units: i32, hom_bases: i32) -> ContextType {
    if repeat_units >= REPETITIVE_REGION_MIN_REPEAT_UNITS {
        if indel_length < SHORT_REPEAT_UNIT_MAX_LENGTH { ContextType::Repeat } else { ContextType::Microhomology }
    } else if hom_bases >= MIN_HOMOLOGOUS_BASES {
        ContextType::Microhomology
    } else if hom_bases >= LONG_INDEL_MIN_HOMOLOGOUS_BASES && indel_length >= LONG_INDEL_MIN_LENGTH {
        ContextType::Microhomology
    } else {
        ContextType::None
    }
}

impl IndelDetails {
    /// IndelDetails.from
    pub fn from(indel: &IndelVariant, ref_genome: &RefGenomeSource) -> Result<IndelDetails> {
        let indel_length = indel.indel_length;
        let mutation_type = if indel.is_deletion { MutationType::Deletion } else { MutationType::Insertion };
        let indel_sequence = indel.indel_sequence();
        let right_flank = indel.flanking_ref_bases_right(ref_genome, RIGHT_FLANK_INDEL_LENGTH_UNITS)?;
        let left_flank = indel.flanking_ref_bases_left(ref_genome, LEFT_FLANK_INDEL_LENGTH_UNITS)?;
        let right_hom = count_homologous_bases(&indel_sequence, &right_flank);
        let left_hom = count_homologous_bases(&reverse(&indel_sequence), &reverse(&left_flank));
        let hom = right_hom.max(left_hom);
        let repeat_units = count_repeat_units(&indel_sequence, &right_flank);
        let context_type = determine_indel_context(indel_length, repeat_units, hom);
        Ok(IndelDetails {
            chromosome: indel.variant.chromosome.clone(),
            position: indel.variant.position,
            ref_bases: indel.variant.ref_bases.clone(),
            alt_bases: indel.variant.alt_bases.clone(),
            indel_sequence,
            indel_length,
            right_flank_bases: right_flank,
            left_flank_bases: left_flank,
            right_hom_bases_count: right_hom,
            left_hom_bases_count: left_hom,
            max_hom_bases_count: hom,
            repeat_units_count: repeat_units,
            mutation_type,
            context_type,
        })
    }

    /// IndelContext.getContextName
    pub fn context_name(&self) -> String {
        let capped = if self.context_type == ContextType::Microhomology {
            self.max_hom_bases_count.min(HOM_BASES_COUNT_CAP)
        } else {
            self.indel_length.min(INDEL_LENGTH_CAP)
        };
        format!("{}{}{}", self.mutation_type.abbreviation(), self.context_type.suffix(), capped)
    }

    pub fn header() -> String {
        [
            "SampleId",
            "Chromosome",
            "Position",
            "RefBases",
            "AltBases",
            "IndelSequence",
            "IndelLength",
            "RightFlankBases",
            "LeftFlankBases",
            "RightHomBasesCount",
            "LeftHomBasesCount",
            "MaxHomBasesCount",
            "RepeatUnitsCount",
            "MutationType",
            "ContextType",
        ]
        .join("\t")
    }

    pub fn line(&self, sample_id: &str) -> String {
        [
            sample_id.to_string(),
            self.chromosome.clone(),
            self.position.to_string(),
            self.ref_bases.clone(),
            self.alt_bases.clone(),
            self.indel_sequence.clone(),
            self.indel_length.to_string(),
            self.right_flank_bases.clone(),
            self.left_flank_bases.clone(),
            self.right_hom_bases_count.to_string(),
            self.left_hom_bases_count.to_string(),
            self.max_hom_bases_count.to_string(),
            self.repeat_units_count.to_string(),
            self.mutation_type.name().to_string(),
            self.context_type.name().to_string(),
        ]
        .join("\t")
    }
}

/// IndelContext.initializeCounts: the 30 bin names in output order
pub fn context_names() -> Vec<String> {
    let mut out = Vec::new();
    for c in ContextType::VALUES {
        for m in MutationType::VALUES {
            let cap = if c == ContextType::Microhomology { HOM_BASES_COUNT_CAP } else { INDEL_LENGTH_CAP };
            for i in 1..=cap {
                out.push(format!("{}{}{}", m.abbreviation(), c.suffix(), i));
            }
        }
    }
    out
}

/// IndelPrep.loadVariants: the INDEL records (VariantContext.isIndel)
pub fn load_indels(variants: &[VcfRecord], log_prefix: &str) -> Vec<IndelVariant> {
    let mut skipped: Vec<(VariantType, i32)> = Vec::new();
    let mut out = Vec::new();
    for v in variants {
        let t = v.variant_type();
        if t != VariantType::Indel {
            match skipped.iter_mut().find(|(k, _)| *k == t) {
                Some(e) => e.1 += 1,
                None => skipped.push((t, 1)),
            }
            continue;
        }
        out.push(IndelVariant::new(SmallVariant::from(v)));
    }
    debug!("{log_prefix}Found {} indels", out.len());
    if !skipped.is_empty() {
        let text: Vec<String> = skipped.iter().map(|(k, n)| format!("{}={n}", type_name(*k))).collect();
        debug!("{log_prefix}Skipped variants: {{{}}}", text.join(", "));
    }
    out
}

/// IndelPrep.countMutationContexts after loading; `write_details` gets the details of every indel
pub fn count_mutation_contexts(
    indels: &[IndelVariant],
    ref_genome: &LazyRefGenome,
    keep_details: bool,
    log_prefix: &str,
    mut write_details: impl FnMut(&[IndelDetails]) -> Result<()>,
) -> Result<Vec<MutContextCount>> {
    let names = context_names();
    let mut counts = vec![0i32; names.len()];
    let mut details = Vec::new();
    for indel in indels {
        let d = IndelDetails::from(indel, ref_genome.get()?)?;
        let name = d.context_name();
        let Some(i) = names.iter().position(|n| *n == name) else {
            bail!("java.lang.NullPointerException: no indel context {name}");
        };
        counts[i] += 1;
        if keep_details {
            details.push(d);
        }
    }
    if keep_details {
        write_details(&details)?;
    }
    debug!("{log_prefix}Completed IndelPrep");
    Ok(names.into_iter().zip(counts).map(|(name, count)| MutContextCount { name, count }).collect())
}
