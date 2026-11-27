//! HgvsProtein: the HGVS protein (p.) string of a coding impact.

use anyhow::Result;

use super::effects::{is_inframe, is_nonsense_or_frameshift, VariantEffect};
use super::hgvs_coding::HGVS_UNKNOWN;
use super::{char_at, ProteinContext, SE_END, SE_START};

pub const HGVS_SPLICE_UNKNOWN: &str = "p.?";
const HGVS_FRAMESHIFT: &str = "fs";
const HGVS_STOP_LOST: &str = "ext*?";
const HGVS_STOP_GAINED: &str = "*";
const HGVS_STOP_TRI_CODE: &str = "Ter";
const HGVS_SYNONYMOUS: &str = "=";

/// AminoAcids.TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER reversed, with the stop as Ter
fn tri_letters(aa: u8) -> &'static str {
    match aa {
        b'A' => "Ala",
        b'C' => "Cys",
        b'D' => "Asp",
        b'E' => "Glu",
        b'F' => "Phe",
        b'G' => "Gly",
        b'H' => "His",
        b'I' => "Ile",
        b'K' => "Lys",
        b'L' => "Leu",
        b'M' => "Met",
        b'N' => "Asn",
        b'P' => "Pro",
        b'Q' => "Gln",
        b'R' => "Arg",
        b'S' => "Ser",
        b'T' => "Thr",
        b'V' => "Val",
        b'W' => "Trp",
        b'Y' => "Tyr",
        b'X' => HGVS_STOP_TRI_CODE,
        _ => HGVS_UNKNOWN,
    }
}

fn tri_letters_str(aas: &str) -> String {
    aas.bytes().map(tri_letters).collect()
}

pub fn report_protein_impact(effect: VariantEffect) -> bool {
    is_inframe(effect) || is_nonsense_or_frameshift(effect) || effect == VariantEffect::Synonymous || effect == VariantEffect::Missense
}

/// generate(proteinContext, effects): from the first effect other than stop-gained, a stop gained marked
pub fn generate(pc: &ProteinContext, effects: &[VariantEffect]) -> Result<String> {
    let mut has_stop_gained = false;
    let mut top_effect: Option<VariantEffect> = None;
    for &e in effects {
        if e == VariantEffect::StopGained {
            has_stop_gained = true;
        } else if top_effect.is_none() {
            top_effect = Some(e);
        }
    }
    if top_effect.is_none() && has_stop_gained {
        // reinstated for the protein string if the effect was removed
        top_effect = Some(VariantEffect::Missense);
    }
    let Some(top_effect) = top_effect else { anyhow::bail!("NullPointerException: no protein effect") };
    let mut hgvs = generate_for(pc, top_effect)?;
    if has_stop_gained {
        if let Some(stop_index) = hgvs.find(HGVS_STOP_TRI_CODE) {
            // replace from Ter onwards
            hgvs = format!("{}{}", &hgvs[..stop_index], HGVS_STOP_GAINED);
        } else if top_effect != VariantEffect::Frameshift {
            hgvs.push_str(HGVS_STOP_GAINED);
        } else if hgvs.ends_with(HGVS_FRAMESHIFT) && char_at(&pc.net_alt_amino_acids, 0)? == b'X' {
            hgvs = format!("{}{}", &hgvs[..hgvs.len() - HGVS_FRAMESHIFT.len()], HGVS_STOP_GAINED);
        }
    }
    Ok(hgvs)
}

/// generate(proteinContext, effect)
pub fn generate_for(pc: &ProteinContext, effect: VariantEffect) -> Result<String> {
    if !pc.valid_ref_codon() {
        return Ok(HGVS_UNKNOWN.to_string());
    }
    if !pc.valid_alt_codon() && effect != VariantEffect::Frameshift {
        // all others should have complete codons
        return Ok(HGVS_UNKNOWN.to_string());
    }
    let mut sb = String::from("p.");
    use VariantEffect::*;
    match effect {
        Synonymous | PhasedSynonymous => form_synonymous(pc, &mut sb)?,
        Missense | PhasedMissense => form_missense(pc, &mut sb)?,
        StartLost => sb.push_str("Met1?"),
        StopLost => form_stop_lost(pc, &mut sb)?,
        InframeDeletion | PhasedInframeDeletion => form_inframe_deletion(pc, &mut sb)?,
        InframeInsertion | PhasedInframeInsertion => form_inframe_insertion(pc, &mut sb)?,
        Frameshift => form_frameshift(pc, &mut sb)?,
        _ => sb.push_str(HGVS_UNKNOWN),
    }
    Ok(sb)
}

/// p.Leu54= for 1 codon or p.Leu54_Arg55= for an MNV or phased synonymous
fn form_synonymous(pc: &ProteinContext, sb: &mut String) -> Result<()> {
    sb.push_str(tri_letters(char_at(&pc.ref_amino_acids, 0)?));
    sb.push_str(&pc.codon_index.to_string());
    if pc.ref_amino_acids.len() > 1 {
        sb.push('_');
        let last = pc.ref_amino_acids.len() as i32 - 1;
        sb.push_str(tri_letters(char_at(&pc.ref_amino_acids, last)?));
        sb.push_str(&(pc.codon_index + last).to_string());
    }
    sb.push_str(HGVS_SYNONYMOUS);
    Ok(())
}

/// from an SNV or MNV: p.Arg64Lys, or p.Ala100_Val101delinsArgTrp
fn form_missense(pc: &ProteinContext, sb: &mut String) -> Result<()> {
    let start = pc.net_codon_index_range[SE_START];
    let end = pc.net_codon_index_range[SE_END];
    let r = &pc.net_ref_amino_acids;
    let a = &pc.net_alt_amino_acids;
    if r.len() == 1 {
        sb.push_str(tri_letters(char_at(r, 0)?));
        sb.push_str(&start.to_string());
        sb.push_str(tri_letters(char_at(a, 0)?));
    } else {
        sb.push_str(tri_letters(char_at(r, 0)?));
        sb.push_str(&start.to_string());
        let last = r.len() as i32 - 1;
        sb.push('_');
        sb.push_str(tri_letters(char_at(r, last)?));
        sb.push_str(&end.to_string());
        sb.push_str("delins");
        for i in 0..=last {
            sb.push_str(tri_letters(char_at(a, i)?));
        }
    }
    Ok(())
}

/// p.Lys2del, p.Gly4_Gln6del, p.Cys28_Lys29delinsTrp
fn form_inframe_deletion(pc: &ProteinContext, sb: &mut String) -> Result<()> {
    if pc.is_phased && pc.net_ref_amino_acids.len() > 1 && pc.net_alt_amino_acids.len() > 1 {
        return form_phased_inframe(pc, sb);
    }
    let start = pc.net_codon_index_range[SE_START];
    let end = pc.net_codon_index_range[SE_END];
    let r = &pc.net_ref_amino_acids;
    let a = &pc.net_alt_amino_acids;
    if r.len() == 1 {
        sb.push_str(tri_letters(char_at(r, 0)?));
        sb.push_str(&start.to_string());
        sb.push_str("del");
    } else {
        sb.push_str(tri_letters(char_at(r, 0)?));
        sb.push_str(&start.to_string());
        sb.push('_');
        sb.push_str(tri_letters(char_at(r, r.len() as i32 - 1)?));
        sb.push_str(&end.to_string());
        sb.push_str("del");
        if !a.is_empty() {
            sb.push_str("ins");
            sb.push_str(tri_letters(char_at(a, 0)?));
        }
    }
    Ok(())
}

/// duplications p.Gln8dup, p.Gly4_Gln6dup; insertions p.A2_C3insR, p.K2_L3insGS, p.C28delinsWV, p.Ser81_Val82ins*
fn form_inframe_insertion(pc: &ProteinContext, sb: &mut String) -> Result<()> {
    if pc.is_phased && pc.net_ref_amino_acids.len() > 1 && pc.net_alt_amino_acids.len() > 1 {
        return form_phased_inframe(pc, sb);
    }
    if pc.is_duplication {
        let a = &pc.net_alt_amino_acids;
        let start = pc.net_codon_index_range[SE_START];
        let end = pc.net_codon_index_range[SE_END];
        if end > start {
            sb.push_str(tri_letters(char_at(a, 0)?));
            sb.push_str(&start.to_string());
            sb.push('_');
            sb.push_str(tri_letters(char_at(a, a.len() as i32 - 1)?));
            sb.push_str(&end.to_string());
        } else {
            sb.push_str(tri_letters(char_at(a, 0)?));
            sb.push_str(&start.to_string());
        }
        sb.push_str("dup");
    } else {
        let r = &pc.ref_amino_acids;
        if pc.net_ref_amino_acids.is_empty() {
            // skip past any redundant ref codons upstream
            let rb = r.as_bytes();
            let ab = pc.alt_amino_acids.as_bytes();
            let limit = (r.len() as i32 - 1).min(ab.len() as i32);
            let mut i: i32 = 0;
            while i < limit {
                if rb[i as usize] != ab[i as usize] {
                    break;
                }
                i += 1;
            }
            i -= 1;
            if i < 0 || i >= r.len() as i32 - 1 {
                sb.push_str(HGVS_UNKNOWN);
                return Ok(());
            }
            let aa_start = pc.codon_index + i;
            let aa_end = aa_start + 1;
            // conservative: an insert only, no deleted amino acids, so quote the range
            sb.push_str(tri_letters(rb[i as usize]));
            sb.push_str(&aa_start.to_string());
            sb.push('_');
            sb.push_str(tri_letters(rb[i as usize + 1]));
            sb.push_str(&aa_end.to_string());
        } else {
            // non-conservative: an amino acid has been deleted
            sb.push_str(tri_letters(char_at(&pc.net_ref_amino_acids, 0)?));
            sb.push_str(&pc.net_codon_index_range[SE_START].to_string());
            sb.push_str("del");
        }
        sb.push_str("ins");
        sb.push_str(&tri_letters_str(&pc.net_alt_amino_acids));
    }
    Ok(())
}

/// several amino acids deleted and changed or inserted
fn form_phased_inframe(pc: &ProteinContext, sb: &mut String) -> Result<()> {
    let r = &pc.net_ref_amino_acids;
    sb.push_str(tri_letters(char_at(r, 0)?));
    sb.push_str(&pc.net_codon_index_range[SE_START].to_string());
    sb.push('_');
    sb.push_str(tri_letters(char_at(r, r.len() as i32 - 1)?));
    sb.push_str(&pc.net_codon_index_range[SE_END].to_string());
    sb.push_str("delins");
    sb.push_str(&tri_letters_str(&pc.net_alt_amino_acids));
    Ok(())
}

/// the first changed (ref) amino acid downstream and its index, otherwise the last ref one
fn form_frameshift(pc: &ProteinContext, sb: &mut String) -> Result<()> {
    let r = pc.ref_amino_acids.as_bytes();
    let a = pc.alt_amino_acids.as_bytes();
    let mut i: i32 = 0;
    while (i as usize) < r.len() {
        if i as usize >= a.len() {
            break;
        }
        if r[i as usize] != a[i as usize] {
            break;
        }
        i += 1;
    }
    if i as usize == r.len() {
        i -= 1;
    }
    sb.push_str(tri_letters(char_at(&pc.ref_amino_acids, i)?));
    sb.push_str(&(pc.codon_index + i).to_string());
    sb.push_str(HGVS_FRAMESHIFT);
    Ok(())
}

/// p.Ter494ext*?, without the new amino acids after the stop
fn form_stop_lost(pc: &ProteinContext, sb: &mut String) -> Result<()> {
    sb.push_str(tri_letters(char_at(&pc.ref_amino_acids, 0)?));
    sb.push_str(&pc.codon_index.to_string());
    let r = pc.ref_amino_acids.as_bytes();
    let a = pc.alt_amino_acids.as_bytes();
    if r.len() == a.len() {
        sb.push_str(tri_letters(char_at(&pc.alt_amino_acids, a.len() as i32 - 1)?));
    } else {
        // a frameshift if the first different amino acid is not the stop
        for i in 0..r.len().min(a.len()) {
            if r[i] != a[i] {
                if a[i] != a[0] {
                    sb.push_str(tri_letters(a[i]));
                }
                sb.push_str(HGVS_FRAMESHIFT);
                return Ok(());
            }
        }
    }
    sb.push_str(HGVS_STOP_LOST);
    Ok(())
}
