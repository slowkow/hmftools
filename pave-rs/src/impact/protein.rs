//! ProteinUtils: the codons and amino acids a coding variant changes.

use anyhow::Result;

use hmf_common::codon::{amino_acid_from_bases, reverse_complement_bases};
use hmf_common::ensembl::{ExonData, TranscriptData};
use hmf_common::regions::position_within;

use super::utils::codon_for_base;
use super::{char_at, coding_end, coding_start, jsub, jsub_from, or_null, CodingContext, ProteinContext, RefGenome, PHASE_1, PHASE_2, SE_END, SE_START};
use crate::variant::VariantData;

const CODON_LENGTH: i32 = 3;

fn aa(bases: &str) -> String {
    amino_acid_from_bases(bases.as_bytes())
}

fn aa_rc(bases: &str) -> String {
    amino_acid_from_bases(&reverse_complement_bases(bases.as_bytes()))
}

fn amino_acids(bases: &str, pos_strand: bool) -> String {
    if pos_strand { aa(bases) } else { aa_rc(bases) }
}

pub fn determine_context(
    variant: &VariantData, cc: &CodingContext, trans: &TranscriptData, ref_genome: &dyn RefGenome,
) -> Result<ProteinContext> {
    let mut pc = ProteinContext::default();
    let pos_strand = trans.pos_strand();

    let Some(exon) = trans.exons.iter().find(|x| x.rank == cc.exon_rank) else {
        crate::error!("var({}) trans({}) invalid coding exonRank({})", variant, trans.trans_name, cc.exon_rank);
        return Ok(pc);
    };

    // MNVs have the coding range exactly over their ref/alt bases (minus non-coding bases), DELs cover the
    // first ref base either side of the DEL, and inserts just the upstream ref base
    pc.codon_index = codon_for_base(cc.coding_base);
    pc.net_codon_index_range = [pc.codon_index, pc.codon_index];

    let ref_coding_bases =
        or_null(ref_genome.base_string(&variant.chromosome, cc.coding_position_range[SE_START], cc.coding_position_range[SE_END])?);
    pc.ref_codon_bases = ref_coding_bases.clone();
    pc.ref_codons_ranges.push([cc.coding_position_range[SE_START], cc.coding_position_range[SE_END]]);

    let upstream_open_codon_bases = if cc.coding_base > 1 { get_open_codon_bases(cc.upstream_phase) } else { 0 };

    if upstream_open_codon_bases > 0 {
        // bases upstream to complete the upstream part of the codon
        let search_up = !pos_strand;
        let upstream =
            or_null(get_extra_bases(trans, ref_genome, &variant.chromosome, exon, upstream_open_codon_bases, search_up, &mut pc.ref_codons_ranges)?);
        if pos_strand {
            pc.ref_codon_bases = upstream + &pc.ref_codon_bases;
        } else {
            pc.ref_codon_bases += &upstream;
        }
    }

    let downstream_mod = pc.ref_codon_bases.len() as i32 % CODON_LENGTH;
    let mut downstream_open_codon_bases = if downstream_mod == 0 { 0 } else { CODON_LENGTH - downstream_mod };
    if variant.is_insert() {
        downstream_open_codon_bases += CODON_LENGTH; // the ref codon past the insert as well
    }

    if downstream_open_codon_bases > 0 {
        let search_up = pos_strand;
        let downstream =
            or_null(get_extra_bases(trans, ref_genome, &variant.chromosome, exon, downstream_open_codon_bases, search_up, &mut pc.ref_codons_ranges)?);
        if pos_strand {
            pc.ref_codon_bases += &downstream;
        } else {
            pc.ref_codon_bases = downstream + &pc.ref_codon_bases;
        }
    }

    // incomplete codons come from a bug or an issue with the transcript definition
    if !pc.valid_ref_codon() {
        return Ok(pc);
    }

    // only the coding bases of the mutation
    let ref_ = cc.coding_ref(variant)?;
    let mut alt = cc.coding_alt(variant)?;
    let var_length = ref_.len() as i32;

    if pos_strand {
        if upstream_open_codon_bases > 0 {
            pc.alt_codon_bases = format!(
                "{}{}{}",
                jsub(&pc.ref_codon_bases, 0, upstream_open_codon_bases)?,
                alt,
                jsub_from(&pc.ref_codon_bases, upstream_open_codon_bases + var_length)?
            );
        } else {
            pc.alt_codon_bases = format!("{}{}", alt, jsub_from(&pc.ref_codon_bases, var_length.min(pc.ref_codon_bases.len() as i32))?);
        }
    } else {
        let codon_base_length = pc.ref_codon_bases.len() as i32;
        // a DEL on the negative strand ending at the exon boundary or coding region has already been
        // truncated by a base
        let upstream_ref_bases = if variant.is_deletion() && cc.coding_position_range[SE_END] == variant.end_position {
            upstream_open_codon_bases + 1
        } else {
            upstream_open_codon_bases
        };
        if variant.is_insert() {
            // inserts have the ref on the upper side
            alt = format!("{}{}", jsub_from(&alt, 1)?, ref_coding_bases);
        } else if variant.is_deletion() {
            alt = if !alt.is_empty() { jsub(&alt, 0, 1)?.to_string() } else { String::new() };
        }
        if upstream_ref_bases > 0 && codon_base_length >= var_length + upstream_ref_bases {
            pc.alt_codon_bases = format!(
                "{}{}{}",
                jsub(&pc.ref_codon_bases, 0, codon_base_length - var_length - upstream_ref_bases)?,
                alt,
                jsub_from(&pc.ref_codon_bases, codon_base_length - upstream_ref_bases)?
            );
        } else {
            pc.alt_codon_bases = format!("{}{}", jsub(&pc.ref_codon_bases, 0, codon_base_length - var_length)?, alt);
        }
    }

    pc.ref_amino_acids = amino_acids(&pc.ref_codon_bases, pos_strand);
    pc.alt_codon_bases_complete = pc.alt_codon_bases.clone();

    // fill in any incomplete alt amino acids from a frameshift
    let downstream_alt_mod = pc.alt_codon_bases.len() as i32 % CODON_LENGTH;
    let downstream_alt_open_codon_bases = if downstream_alt_mod == 0 { 0 } else { CODON_LENGTH - downstream_alt_mod };

    if downstream_alt_open_codon_bases > 0 {
        let downstream_start_pos = if pos_strand { pc.ref_coding_base_end() } else { pc.ref_coding_base_start() };
        let search_up = pos_strand;
        let downstream = or_null(get_extra_bases_from(
            trans,
            ref_genome,
            &variant.chromosome,
            exon,
            downstream_start_pos,
            downstream_alt_open_codon_bases,
            search_up,
        )?);
        if pos_strand {
            pc.alt_codon_bases_complete += &downstream;
        } else {
            pc.alt_codon_bases_complete = downstream + &pc.alt_codon_bases_complete;
        }
    }

    if !pc.valid_alt_codon() {
        pc.net_ref_amino_acids = pc.ref_amino_acids.clone();
        pc.net_codon_index_range[SE_START] = pc.codon_index;
        pc.net_codon_index_range[SE_END] = pc.codon_index + pc.ref_amino_acids.len() as i32 - 1;
        return Ok(pc);
    }

    pc.alt_amino_acids = amino_acids(&pc.alt_codon_bases_complete, pos_strand);

    // extend the up and downstream codons of INDELs to find the net differences in amino acids
    if variant.is_indel() {
        extend_indel_codon_bases(variant, cc, &mut pc, trans, exon, ref_genome, downstream_alt_open_codon_bases)?;
        extend_ref_amino_acids(variant, &mut pc, trans, exon, ref_genome)?;
    }

    trim_variant_amino_acids(variant, &mut pc, pos_strand);

    check_insert_duplication(variant, cc, &mut pc)?;

    Ok(pc)
}

fn extend_indel_codon_bases(
    variant: &VariantData, cc: &CodingContext, pc: &mut ProteinContext, trans: &TranscriptData, exon: &ExonData,
    ref_genome: &dyn RefGenome, downstream_alt_open_codon_bases: i32,
) -> Result<()> {
    let pos_strand = trans.pos_strand();

    // for indels within a codon (phase 1 or 2), take an extra upstream codon
    if variant.is_insert() && (cc.upstream_phase == PHASE_1 || cc.upstream_phase == PHASE_2) && pc.codon_index > 1 {
        // only if the first 1 or 2 ref amino acids match the last 1 or 2 of the alt
        let first_matches_last = char_at(&pc.ref_amino_acids, 0)? == char_at(&pc.alt_amino_acids, pc.alt_amino_acids.len() as i32 - 1)?;
        let first_seq_matches_last = pc.ref_amino_acids.len() >= 2
            && pc.alt_amino_acids.len() >= 2
            && pc.ref_amino_acids[..2] == pc.alt_amino_acids[pc.alt_amino_acids.len() - 2..];

        if first_matches_last || first_seq_matches_last {
            let search_up = !pos_strand;
            let upstream = get_extra_bases(trans, ref_genome, &variant.chromosome, exon, CODON_LENGTH, search_up, &mut pc.ref_codons_ranges)?;
            if let Some(up) = upstream {
                pc.codon_index -= 1;
                pc.extra_upstream_codon = true;
                if pos_strand {
                    pc.ref_codon_bases = up.clone() + &pc.ref_codon_bases;
                    pc.alt_codon_bases = up.clone() + &pc.alt_codon_bases;
                    pc.alt_codon_bases_complete = up + &pc.alt_codon_bases_complete;
                } else {
                    pc.ref_codon_bases += &up;
                    pc.alt_codon_bases += &up;
                    pc.alt_codon_bases_complete += &up;
                }
                pc.ref_amino_acids = amino_acids(&pc.ref_codon_bases, pos_strand);
                pc.alt_amino_acids = amino_acids(&pc.alt_codon_bases_complete, pos_strand);
            }
        }
    }

    if pc.ref_amino_acids != pc.alt_amino_acids || downstream_alt_open_codon_bases == 0 {
        return Ok(());
    }

    // from the last ref codon base take 3 more each time
    let search_up = pos_strand;
    let current_pos = if pos_strand { pc.ref_coding_base_end() } else { pc.ref_coding_base_start() };
    let mut extra_bases = CODON_LENGTH + downstream_alt_open_codon_bases;

    loop {
        let Some(downstream) = get_extra_bases_from(trans, ref_genome, &variant.chromosome, exon, current_pos, extra_bases, search_up)? else {
            break;
        };
        let len = downstream.len() as i32;
        let coding_region_ended = len < extra_bases || len < downstream_alt_open_codon_bases;
        if coding_region_ended {
            // the end of the coding region is reached
            crate::debug!(
                "var({}) invalid extraBases({}) downstreamBases({}) downstreamAltOpenCodonBases({})",
                variant,
                extra_bases,
                downstream,
                downstream_alt_open_codon_bases
            );
        }
        if pos_strand {
            // trim off the extra bases the frame-shifted alt needs
            let ref_downstream = jsub(&downstream, 0, len - downstream_alt_open_codon_bases)?;
            pc.ref_codon_bases_extended = format!("{}{}", pc.ref_codon_bases, ref_downstream);
            pc.alt_codon_bases_complete = format!("{}{}", pc.alt_codon_bases, downstream);
            pc.ref_amino_acids = aa(&pc.ref_codon_bases_extended);
            pc.alt_amino_acids = aa(&pc.alt_codon_bases_complete);
        } else {
            let ref_downstream = jsub_from(&downstream, downstream_alt_open_codon_bases)?;
            pc.ref_codon_bases_extended = format!("{}{}", ref_downstream, pc.ref_codon_bases);
            pc.alt_codon_bases_complete = format!("{}{}", downstream, pc.alt_codon_bases);
            pc.ref_amino_acids = aa_rc(&pc.ref_codon_bases_extended);
            pc.alt_amino_acids = aa_rc(&pc.alt_codon_bases_complete);
        }
        if pc.ref_amino_acids != pc.alt_amino_acids {
            break;
        }
        if coding_region_ended {
            break;
        }
        extra_bases += CODON_LENGTH;
    }
    Ok(())
}

fn extend_ref_amino_acids(
    variant: &VariantData, pc: &mut ProteinContext, trans: &TranscriptData, exon: &ExonData, ref_genome: &dyn RefGenome,
) -> Result<()> {
    if !trans.pos_strand() || !variant.is_insert() {
        return Ok(());
    }

    let current_pos = pc.ref_coding_base_end();
    let mut extra_bases = CODON_LENGTH;

    // if the ref and alt amino acids match, continue as long as the next ref ones match the additional alt
    // ones (a duplication), eg MAPPPPQ -> MAPPPP P Q, or LGHEX -> LGHE GHE X
    let r = pc.ref_amino_acids.as_bytes();
    let a = pc.alt_amino_acids.as_bytes();
    let mut alt_aa_index: i32 = -1;
    for i in 0..r.len().min(a.len()) {
        if r[i] == a[i] {
            alt_aa_index = i as i32;
        } else {
            break;
        }
    }
    if alt_aa_index < 0 {
        return Ok(());
    }

    // search for the different amino acid in the ref from this point
    alt_aa_index += 1;
    let alt_len = a.len() as i32;
    let mut comparison_amino_acid = if alt_aa_index >= alt_len { a[a.len() - 1] } else { a[alt_aa_index as usize] };

    let mut extended_ref_bases = String::new();

    loop {
        let downstream = get_extra_bases_from(trans, ref_genome, &variant.chromosome, exon, current_pos, extra_bases, true)?;
        let Some(downstream) = downstream else { break };
        if (downstream.len() as i32) < CODON_LENGTH {
            break;
        }
        let extended_amino_acids = aa(&downstream);
        let next_ref_amino_acid = *extended_amino_acids.as_bytes().last().unwrap();
        if next_ref_amino_acid != comparison_amino_acid {
            break;
        }
        if alt_aa_index < alt_len - 1 {
            alt_aa_index += 1;
            comparison_amino_acid = a[alt_aa_index as usize];
        }
        let coding_region_ended = (downstream.len() as i32) < extra_bases;
        extended_ref_bases = downstream;
        if coding_region_ended {
            break;
        }
        extra_bases += CODON_LENGTH;
    }

    let initial_alt_codon_bases_complete = pc.alt_codon_bases_complete.clone();

    if !extended_ref_bases.is_empty() {
        pc.ref_codon_bases_extended = format!("{}{}", pc.ref_codon_bases, extended_ref_bases);
        pc.alt_codon_bases_complete = format!("{}{}", initial_alt_codon_bases_complete, extended_ref_bases);
        pc.ref_amino_acids = aa(&pc.ref_codon_bases_extended);
        pc.alt_amino_acids = aa(&pc.alt_codon_bases_complete);
    }

    // ensure there is an extra ref amino acid if the ref ones match the alt ones
    if pc.ref_amino_acids.len() < pc.alt_amino_acids.len() && pc.alt_amino_acids.starts_with(pc.ref_amino_acids.as_str()) {
        // extend with one more ref amino acid
        if let Some(downstream) = get_extra_bases_from(trans, ref_genome, &variant.chromosome, exon, current_pos, extra_bases, true)? {
            pc.ref_codon_bases_extended = format!("{}{}", pc.ref_codon_bases, downstream);
            pc.alt_codon_bases_complete = format!("{}{}", initial_alt_codon_bases_complete, downstream);
            pc.ref_amino_acids = aa(&pc.ref_codon_bases_extended);
            pc.alt_amino_acids = aa(&pc.alt_codon_bases_complete);
        }
    }
    Ok(())
}

/// Trims matching amino acids from the start and end to align any differences to the most 3' end in
/// protein space
fn trim_variant_amino_acids(variant: &VariantData, pc: &mut ProteinContext, pos_strand: bool) {
    if pc.ref_amino_acids == pc.alt_amino_acids {
        return; // synonymous
    }
    let can_trim_start = variant.is_base_change() || variant.is_deletion() || variant.is_insert();
    let repeat_start_removal = can_trim_start && variant.is_indel();
    trim_amino_acids(pc, pos_strand, can_trim_start, repeat_start_removal, true);
}

/// ProteinUtils.trimAminoAcids
pub fn trim_amino_acids(pc: &mut ProteinContext, _pos_strand: bool, can_trim_start: bool, repeat_start_removal: bool, repeat_end_removal: bool) {
    let mut aa_index = pc.codon_index;
    let mut net_ref: &str = &pc.ref_amino_acids;
    let mut net_alt: &str = &pc.alt_amino_acids;

    // strip off the initial amino acid if it doesn't relate to any altered base, as for a DEL whose ref is
    // the last base of a codon
    if can_trim_start {
        loop {
            // strip a matching amino acid from the start if the ref has more than 1
            if net_ref.len() > 1 && !net_alt.is_empty() && net_ref.as_bytes()[0] == net_alt.as_bytes()[0] {
                net_ref = &net_ref[1..];
                net_alt = &net_alt[1..];
                aa_index += 1;
                if !repeat_start_removal {
                    break;
                }
            } else {
                break;
            }
        }
    }

    // and strip matching amino acids from the end
    loop {
        if !net_ref.is_empty() && !net_alt.is_empty() && net_ref.as_bytes()[net_ref.len() - 1] == net_alt.as_bytes()[net_alt.len() - 1] {
            net_ref = &net_ref[..net_ref.len() - 1];
            net_alt = &net_alt[..net_alt.len() - 1];
            if !repeat_end_removal {
                break;
            }
        } else {
            break;
        }
    }

    let (net_ref, net_alt) = (net_ref.to_string(), net_alt.to_string());
    if net_ref.is_empty() {
        aa_index -= 1;
    }
    pc.net_codon_index_range = [aa_index, aa_index + (net_ref.len() as i32 - 1).max(0)];
    pc.net_ref_amino_acids = net_ref;
    pc.net_alt_amino_acids = net_alt;
}

/// ProteinUtils.checkInsertDuplication
pub fn check_insert_duplication(variant: &VariantData, cc: &CodingContext, pc: &mut ProteinContext) -> Result<()> {
    if !variant.is_insert() || cc.is_frameshift {
        return Ok(());
    }
    if pc.net_alt_amino_acids.is_empty() || pc.ref_amino_acids.len() < 2 {
        return Ok(());
    }

    // the first inserted amino acid matches a ref one
    let first_alt = pc.net_alt_amino_acids.as_bytes()[0];
    if !pc.ref_amino_acids.as_bytes().contains(&first_alt) {
        return Ok(());
    }

    // a duplication of amino acids means 1 extra copy: p.Gln8dup, p.Gly4_Gln6dup, c.1014_1019dupTGCTGC p.Ala339_Ala340dup
    if pc.net_alt_amino_acids.len() == 1 && pc.net_ref_amino_acids.is_empty() {
        if first_alt == pc.ref_amino_acids.as_bytes()[0] {
            pc.is_duplication = true;
            pc.net_codon_index_range = [pc.codon_index, pc.codon_index];
            return Ok(());
        } else if first_alt == pc.ref_amino_acids.as_bytes()[1] {
            pc.is_duplication = true;
            pc.net_codon_index_range[SE_END] = pc.net_codon_index_range[SE_START];
            return Ok(());
        }
    }

    // no partial duplications, so find the index in the ref
    let Some(dup_index) = pc.ref_amino_acids.find(pc.net_alt_amino_acids.as_str()) else { return Ok(()) };
    pc.is_duplication = true;
    pc.net_codon_index_range[SE_START] = pc.codon_index + dup_index as i32;
    pc.net_codon_index_range[SE_END] = pc.net_codon_index_range[SE_START] + pc.net_alt_amino_acids.len() as i32 - 1;
    Ok(())
}

/// the number of bases upstream to complete a codon
pub fn get_open_codon_bases(phase: i32) -> i32 {
    if phase == PHASE_1 {
        0
    } else if phase == PHASE_2 {
        1
    } else {
        2
    }
}

/// getExtraBases(transData, refGenome, chromosome, currentExon, startPos, requiredBases, searchUp)
pub fn get_extra_bases_from(
    trans: &TranscriptData, ref_genome: &dyn RefGenome, chromosome: &str, exon: &ExonData, start_pos: i32, required_bases: i32,
    search_up: bool,
) -> Result<Option<String>> {
    let mut ranges = vec![[start_pos, start_pos]];
    get_extra_bases(trans, ref_genome, chromosome, exon, required_bases, search_up, &mut ranges)
}

/// getExtraBases: `required_bases` further coding bases from the end (search up) or start of the coding
/// base ranges, which are extended; None (Java's null) when the start is in no exon
pub fn get_extra_bases(
    trans: &TranscriptData, ref_genome: &dyn RefGenome, chromosome: &str, ref_exon: &ExonData, mut required_bases: i32,
    search_up: bool, ranges: &mut Vec<[i32; 2]>,
) -> Result<Option<String>> {
    if ranges.is_empty() {
        return Ok(None);
    }
    let current_index = if search_up { ranges.len() - 1 } else { 0 };
    let start_pos = if search_up { ranges[current_index][SE_END] } else { ranges[current_index][SE_START] };

    let current_exon = if position_within(start_pos, ref_exon.start, ref_exon.end) {
        Some(ref_exon)
    } else {
        trans.exons.iter().find(|x| position_within(start_pos, x.start, x.end))
    };
    let Some(current_exon) = current_exon else { return Ok(None) };

    // Java's String concatenation: a null String appends "null"
    let mut extra_bases: Option<String> = Some(String::new());

    if search_up {
        let current_exon_bases = required_bases.min(current_exon.end.min(coding_end(trans)?) - start_pos);
        if current_exon_bases > 0 {
            let upper_pos = start_pos + 1 + (current_exon_bases - 1);
            ranges[current_index][SE_END] = upper_pos;
            extra_bases = ref_genome.base_string(chromosome, start_pos + 1, upper_pos)?;
            required_bases -= current_exon_bases;
        }
        if required_bases > 0 {
            let next_rank = if trans.pos_strand() { current_exon.rank + 1 } else { current_exon.rank - 1 };
            let next_exon = trans.exons.iter().find(|x| x.rank == next_rank);
            let Some(next_exon) = next_exon.filter(|ne| ne.start <= coding_end(trans).unwrap_or(i32::MAX)) else { return Ok(extra_bases) };
            let lower = next_exon.start;
            let upper = next_exon.start + (required_bases - 1);
            ranges.push([lower, upper]);
            let next_bases = or_null(ref_genome.base_string(chromosome, lower, upper)?);
            extra_bases = Some(or_null(extra_bases) + &next_bases);
        }
    } else {
        // search lower positions
        let current_exon_bases = required_bases.min(start_pos - current_exon.start.max(coding_start(trans)?));
        if current_exon_bases > 0 {
            let lower_pos = start_pos - current_exon_bases;
            ranges[current_index][SE_START] = lower_pos;
            extra_bases = ref_genome.base_string(chromosome, lower_pos, start_pos - 1)?;
            required_bases -= current_exon_bases;
        }
        if required_bases > 0 {
            let next_rank = if trans.pos_strand() { current_exon.rank - 1 } else { current_exon.rank + 1 };
            let next_exon = trans.exons.iter().find(|x| x.rank == next_rank);
            let Some(next_exon) = next_exon.filter(|ne| ne.end >= coding_start(trans).unwrap_or(i32::MIN)) else { return Ok(extra_bases) };
            let lower = next_exon.end - (required_bases - 1);
            let upper = next_exon.end;
            ranges.insert(0, [lower, upper]);
            let next_bases = or_null(ref_genome.base_string(chromosome, lower, upper)?);
            extra_bases = Some(next_bases + &or_null(extra_bases));
        }
    }
    Ok(extra_bases)
}
