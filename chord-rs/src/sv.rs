//! SvPrep, StructuralVariant, SvContext and SvDetails: structural variants (one record per breakend pair)
//! counted by type and length: DEL, DUP and INV in six length bins each, and translocations (BND).

use anyhow::{bail, Result};
use std::collections::HashSet;

use hmf_common::debug;
use hmf_common::sv::{StructuralVariantType, VariantAltInsertCoords, MATE_ID};
use hmf_common::variant_context::{VariantType, VcfRecord};

use crate::compat;
use crate::prep::MutContextCount;
use crate::vcf_file::type_name;

pub const SV_DETAILS_FILE_SUFFIX: &str = ".chord.sv.details.tsv";

/// VariantContext.MAX_ALLELE_SIZE_FOR_NON_SV
const MAX_ALLELE_SIZE_FOR_NON_SV: i32 = 150;

/// SvContext.SV_LENGTH_INTERVALS
const SV_LENGTH_INTERVALS: [i64; 7] = [0, 1_000, 10_000, 100_000, 1_000_000, 10_000_000, i32::MAX as i64];
/// DecimalFormat("0E00") of the intervals, with a lower-case 'e' ("Inf" for Integer.MAX_VALUE)
const SV_LENGTH_NAMES: [&str; 7] = ["0e00", "1e03", "1e04", "1e05", "1e06", "1e07", "Inf"];
const TRANSLOCATION: &str = "TRA";

/// StructuralVariant: one breakend and the type and length it gives with its mate
#[derive(Clone, Debug)]
pub struct StructuralVariant {
    pub id: String,
    pub ref_chromosome: String,
    pub ref_position: i32,
    pub ref_orientation: Option<i8>,
    pub ref_string: String,
    pub alt_string: String,
    pub alt_chromosome: Option<String>,
    pub alt_position: i32,
    pub alt_orientation: Option<i8>,
    pub sv_type: StructuralVariantType,
    pub length: Option<i32>,
}

impl StructuralVariant {
    pub fn new(id: &str, ref_chromosome: &str, ref_position: i32, ref_sequence: &str, alt_string: &str) -> Result<StructuralVariant> {
        let parsed = VariantAltInsertCoords::from_ref_alt_split(alt_string, ref_sequence, !compat::java("CH8"))?;
        let alt_chromosome = if parsed.other_chromosome.is_empty() { None } else { Some(parsed.other_chromosome.clone()) };
        let mut sv = StructuralVariant {
            id: id.to_string(),
            ref_chromosome: ref_chromosome.to_string(),
            ref_position,
            ref_orientation: parsed.orient,
            ref_string: ref_sequence.to_string(),
            alt_string: alt_string.to_string(),
            alt_chromosome,
            alt_position: parsed.other_position,
            alt_orientation: parsed.other_orient,
            sv_type: StructuralVariantType::SGL,
            length: None,
        };
        sv.sv_type = sv.get_type();
        sv.length = match sv.sv_type {
            StructuralVariantType::DEL | StructuralVariantType::DUP | StructuralVariantType::INV | StructuralVariantType::INS => {
                Some((sv.alt_position as i64 - sv.ref_position as i64).abs() as i32)
            }
            _ => None,
        };
        Ok(sv)
    }

    pub fn from(rec: &VcfRecord) -> Result<StructuralVariant> {
        StructuralVariant::new(&rec.id, &rec.chrom, rec.pos, &rec.ref_allele, rec.alt())
    }

    /// StructuralVariant.getType
    fn get_type(&self) -> StructuralVariantType {
        use StructuralVariantType::*;
        let Some(alt_chromosome) = &self.alt_chromosome else { return SGL };
        if &self.ref_chromosome != alt_chromosome {
            return BND;
        }
        let pos_diff = (self.alt_position as i64 - self.ref_position as i64).abs();
        // CHORD was trained with this definition of an insertion (ESVEE also asks for inserted bases)
        if pos_diff == 1 {
            return INS;
        }
        if self.ref_orientation == self.alt_orientation {
            return INV;
        }
        let forward = self.ref_orientation == Some(1);
        let reverse = self.ref_orientation == Some(-1);
        if self.ref_position < self.alt_position {
            if forward { DEL } else { DUP }
        } else if reverse {
            DEL
        } else {
            DUP
        }
    }

    /// SvContext.getContextName
    pub fn context_name(&self) -> String {
        if self.sv_type == StructuralVariantType::BND {
            return TRANSLOCATION.to_string();
        }
        let length = self.length.unwrap_or(0) as i64;
        let i = (1..SV_LENGTH_INTERVALS.len()).find(|&i| length < SV_LENGTH_INTERVALS[i]).unwrap_or(SV_LENGTH_INTERVALS.len() - 1);
        format!("{}_{}_{}_bp", self.sv_type.name(), SV_LENGTH_NAMES[i - 1], SV_LENGTH_NAMES[i])
    }

    pub fn header() -> String {
        ["Id", "RefChromosome", "RefPosition", "RefOrientation", "RefString", "AltChromosome", "AltPosition", "AltOrientation", "AltString", "Type", "Length"]
            .join("\t")
    }

    /// SvDetails.writeLine (only DEL, DUP, INV and BND reach it, which have orientations)
    pub fn line(&self) -> String {
        let orient = |o: Option<i8>| o.map(|v| v.to_string()).unwrap_or_else(|| "null".to_string());
        [
            self.id.clone(),
            self.ref_chromosome.clone(),
            self.ref_position.to_string(),
            orient(self.ref_orientation),
            self.ref_string.clone(),
            self.alt_chromosome.clone().unwrap_or_else(|| "null".to_string()),
            self.alt_position.to_string(),
            orient(self.alt_orientation),
            self.alt_string.clone(),
            self.sv_type.name().to_string(),
            self.length.map(|l| l.to_string()).unwrap_or_else(|| "null".to_string()),
        ]
        .join("\t")
    }
}

/// SvContext.initializeCounts: DEL, DUP and INV by length, then TRA
pub fn context_names() -> Vec<String> {
    let mut out = Vec::new();
    for t in ["DEL", "DUP", "INV"] {
        for i in 0..SV_LENGTH_INTERVALS.len() - 1 {
            out.push(format!("{t}_{}_{}_bp", SV_LENGTH_NAMES[i], SV_LENGTH_NAMES[i + 1]));
        }
    }
    out.push(TRANSLOCATION.to_string());
    out
}

/// VariantContext.isSymbolicOrSV: a symbolic record, or an INDEL inserting more than 150 bases
fn is_symbolic_or_sv(rec: &VcfRecord) -> bool {
    match rec.variant_type() {
        VariantType::Symbolic => true,
        VariantType::Indel => rec.alts.iter().any(|a| a.len() as i32 - rec.ref_allele.len() as i32 > MAX_ALLELE_SIZE_FOR_NON_SV),
        _ => false,
    }
}

/// SvPrep.loadVariants: the first breakend of each pair (by ID and MATEID), DEL, DUP, INV and BND only
pub fn load_svs(variants: &[VcfRecord], log_prefix: &str) -> Result<Vec<StructuralVariant>> {
    let mut ids: HashSet<String> = HashSet::new();
    let mut out = Vec::new();
    for rec in variants {
        if !is_symbolic_or_sv(rec) {
            bail!(
                "java.lang.IllegalStateException: Unexpected non-structural variant: variantType({}) variant({}:{} {}>{})",
                type_name(rec.variant_type()),
                rec.chrom,
                rec.pos,
                rec.ref_allele,
                rec.alts.join(",")
            );
        }
        let id = rec.id.clone();
        // for single breakends, the id stands for the mate id
        let mate_id = rec.info.as_string(MATE_ID, &id);
        if !ids.contains(&id) {
            out.push(StructuralVariant::from(rec)?);
        }
        ids.insert(id);
        ids.insert(mate_id);
    }
    use StructuralVariantType::*;
    let out: Vec<StructuralVariant> = out.into_iter().filter(|v| matches!(v.sv_type, DEL | DUP | INV | BND)).collect();
    debug!("{log_prefix}Found {} SVs", out.len());
    Ok(out)
}

/// SvPrep.countMutationContexts after loading
pub fn count_mutation_contexts(svs: &[StructuralVariant], log_prefix: &str) -> Result<Vec<MutContextCount>> {
    let names = context_names();
    let mut counts = vec![0i32; names.len()];
    for sv in svs {
        let name = sv.context_name();
        let Some(i) = names.iter().position(|n| *n == name) else {
            bail!("java.lang.NullPointerException: no SV context {name}");
        };
        counts[i] += 1;
    }
    debug!("{log_prefix}Completed SvPrep");
    Ok(names.into_iter().zip(counts).map(|(name, count)| MutContextCount { name, count }).collect())
}
