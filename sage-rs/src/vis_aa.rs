//! The visualiser's amino acid track: ports of `vis/AminoAcidEvent` and `vis/AminoAcidUtil`.
//! Errors stand for the exceptions Java throws; the caller then writes no file for the variant, as
//! `VisFileBuilder.writeToHtmlFile` does.

use crate::types::*;
use crate::vis_svg::{GeneRegionKind, GeneRegionViewModel};
use std::collections::VecDeque;

pub const REF_BUFFER_SIZE: i32 = 100;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AminoAcidEvent {
    Substitution { pos: i32, alt: u8 },
    Ins { pos: i32, length: i32 },
    Del { start: i32, end: i32 },
    Ext { pos: i32 },
    StartLost { pos: i32 },
}

impl AminoAcidEvent {
    pub fn start_pos(&self) -> i32 {
        match *self {
            AminoAcidEvent::Substitution { pos, .. }
            | AminoAcidEvent::Ins { pos, .. }
            | AminoAcidEvent::Ext { pos }
            | AminoAcidEvent::StartLost { pos } => pos,
            AminoAcidEvent::Del { start, .. } => start,
        }
    }
    fn end_pos(&self) -> i32 {
        match *self {
            AminoAcidEvent::Del { end, .. } => end,
            _ => self.start_pos(),
        }
    }
    fn priority(&self) -> i32 {
        match self {
            AminoAcidEvent::Substitution { .. } => 4,
            AminoAcidEvent::Ins { .. } => 5,
            AminoAcidEvent::Del { .. } => 3,
            AminoAcidEvent::Ext { .. } => 2,
            AminoAcidEvent::StartLost { .. } => 1,
        }
    }
}

/// AminoAcids.TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER
fn tri_to_single(tri: &str) -> Option<u8> {
    Some(match tri {
        "Ala" => b'A',
        "Cys" => b'C',
        "Asp" => b'D',
        "Glu" => b'E',
        "Phe" => b'F',
        "Gly" => b'G',
        "His" => b'H',
        "Ile" => b'I',
        "Lys" => b'K',
        "Leu" => b'L',
        "Met" => b'M',
        "Asn" => b'N',
        "Pro" => b'P',
        "Gln" => b'Q',
        "Arg" => b'R',
        "Ser" => b'S',
        "Thr" => b'T',
        "Val" => b'V',
        "Trp" => b'W',
        "Tyr" => b'Y',
        _ => return None,
    })
}

fn is_tri(s: &[u8]) -> bool {
    s.len() >= 3 && s[0].is_ascii_uppercase() && s[1].is_ascii_lowercase() && s[2].is_ascii_lowercase()
}

/// the pattern `^([A-Z][a-z][a-z])(\d+)(|[^\d].*)$`: (amino acid, position, length matched)
fn leading_ref(s: &str) -> Option<(&str, i32, usize)> {
    let b = s.as_bytes();
    if !is_tri(b) {
        return None;
    }
    let digits = b[3..].iter().take_while(|c| c.is_ascii_digit()).count();
    if digits == 0 {
        return None;
    }
    Some((&s[..3], s[3..3 + digits].parse().ok()?, 3 + digits))
}

/// AminoAcidEvent.parse
pub fn parse_amino_acid_events(hgvs: &str) -> Result<Vec<AminoAcidEvent>, String> {
    use AminoAcidEvent::*;
    // S38, Java: a stop written as `Ter` in a substitution (p.Gly12Ter) is looked up in the amino acid
    // map and is a null pointer exception; an unknown three-letter code in a `delins` adds the four
    // letters of "null" as amino acids. Corrected: `Ter` is the stop, an unknown code is an error.
    let java_codes = crate::compat::java("S38");
    let bad = || format!("hgvs({hgvs}) does not follow expected format");
    if !hgvs.is_ascii() || hgvs.len() < 2 || &hgvs[..2] != "p." {
        return Err(bad());
    }
    let mut working = &hgvs[2..];
    let (tri1, pos1, len) = leading_ref(working).ok_or_else(bad)?;
    working = &working[len..];
    let mut ref2: Option<i32> = None;
    if working.starts_with('_') {
        if tri1 == "Ter" {
            return Err(bad());
        }
        working = &working[1..];
        let (tri2, pos2, len) = leading_ref(working).ok_or_else(bad)?;
        working = &working[len..];
        if tri2 == "Ter" || pos2 <= pos1 {
            return Err(bad());
        }
        ref2 = Some(pos2);
    }
    if working == "=" {
        return Ok(Vec::new());
    }
    if working == "fs" {
        if ref2.is_some() {
            return Err(bad());
        }
        return Ok(vec![Ext { pos: pos1 }]);
    }
    if tri1 == "Ter" {
        if ref2.is_some() {
            return Err(bad());
        }
        if !(is_tri(working.as_bytes()) && &working[3..] == "ext*?") || &working[..3] == "Ter" {
            return Err(bad());
        }
        return Ok(vec![Ext { pos: pos1 }]);
    }
    if working == "ext*?" {
        if ref2.is_some() {
            return Err(bad());
        }
        return Ok(vec![Ext { pos: pos1 }]);
    }
    if working == "*" || (working.len() == 3 && is_tri(working.as_bytes())) {
        if ref2.is_some() {
            return Err(bad());
        }
        let alt = if working == "*" || (!java_codes && working == "Ter") { b'*' } else { tri_to_single(working).ok_or_else(bad)? };
        return Ok(vec![Substitution { pos: pos1, alt }]);
    }
    if working == "?" {
        if ref2.is_some() {
            return Err(bad());
        }
        return Ok(vec![StartLost { pos: pos1 }]);
    }
    if working == "del" {
        return Ok(vec![Del { start: pos1, end: ref2.unwrap_or(pos1) }]);
    }
    if working == "dup" {
        let end = ref2.unwrap_or(pos1);
        return Ok(vec![Ins { pos: end, length: end - pos1 + 1 }]);
    }
    // Java takes subSequence(0, 3) and then (0, 6), which throw when the text is shorter
    if working.len() < 3 {
        return Err(bad());
    }
    let is_insert = &working[..3] == "ins";
    if !is_insert && (working.len() < 6 || &working[..6] != "delins") {
        return Err(bad());
    }
    working = if is_insert { &working[3..] } else { &working[6..] };
    let ref2 = match ref2 {
        Some(p) => p,
        None if !is_insert => pos1,
        None => return Err(bad()),
    };
    if is_insert && ref2 != pos1 + 1 {
        return Err(bad());
    }
    let mut alt_acids: Vec<u8> = Vec::new();
    while !working.is_empty() {
        if is_tri(working.as_bytes()) {
            let tri = &working[..3];
            if tri == "Ter" {
                return Err(bad());
            }
            match tri_to_single(tri) {
                Some(c) => alt_acids.push(c),
                None if java_codes => alt_acids.extend_from_slice(b"null"),
                None => return Err(bad()),
            }
            working = &working[3..];
            continue;
        }
        if working == "*" {
            alt_acids.push(b'*');
            break;
        }
        return Err(bad());
    }
    let alt_len = alt_acids.len() as i32;
    if is_insert {
        return Ok(vec![Ins { pos: pos1, length: alt_len }]);
    }
    let mut events = Vec::new();
    let gap_size = ref2 - pos1 + 1;
    if alt_len > gap_size {
        let ins_len = alt_len - gap_size;
        events.push(Ins { pos: pos1 - 1, length: ins_len });
        for i in pos1..=ref2 {
            events.push(Substitution { pos: i, alt: alt_acids[(i - pos1 + ins_len) as usize] });
        }
        return Ok(events);
    }
    let del_len = gap_size - alt_len;
    if del_len > 0 {
        events.push(Del { start: pos1, end: pos1 + del_len - 1 });
    }
    for (i, alt) in alt_acids.iter().enumerate() {
        events.push(Substitution { pos: pos1 + i as i32 + del_len, alt: *alt });
    }
    Ok(events)
}

/// The transcript a reportable variant's impact names: TranscriptData with its TranscriptAminoAcids
#[derive(Clone, Debug)]
pub struct Transcript {
    pub gene_id: String,
    pub gene_name: String,
    pub trans_name: String,
    pub pos_strand: bool,
    pub coding_start: Option<i32>,
    pub coding_end: Option<i32>,
    /// (start, end, start phase), in file order
    pub exons: Vec<(i32, i32, i32)>,
    pub amino_acids: Vec<u8>,
}

const PHASE_NONE: i32 = -1;

/// TranscriptUtils.calcExonicCodingPhase
fn calc_exonic_coding_phase(exon: (i32, i32, i32), coding_start: i32, coding_end: i32, pos_strand: bool, position: i32) -> i32 {
    if !position_within(position, coding_start, coding_end) {
        return PHASE_NONE;
    }
    let (exon_start, exon_end, phase_start) = exon;
    let coding_bases =
        if pos_strand { position - exon_start.max(coding_start) + 1 } else { exon_end.min(coding_end) - position + 1 };
    if (pos_strand && coding_start == exon_start) || (!pos_strand && coding_end == exon_end) {
        let start_phase = if phase_start == PHASE_NONE { 0 } else { phase_start };
        (start_phase + coding_bases) % 3
    } else if (pos_strand && coding_start < exon_start) || (!pos_strand && coding_end > exon_end) {
        (phase_start + coding_bases) % 3
    } else {
        coding_bases % 3
    }
}

/// AminoAcidUtil.getGeneRegions: in transcript order; the amino acid letters are filled in later
pub fn gene_regions(t: &Transcript) -> Result<Vec<GeneRegionViewModel>, String> {
    let pos_strand = t.pos_strand;
    let (Some(coding_start), Some(coding_end)) = (t.coding_start, t.coding_end) else {
        return Err(format!("transcript({}) is not coding", t.trans_name));
    };
    let mut regions = Vec::new();
    let mut prev_exon: Option<(i32, i32, i32)> = None;
    let pos_inc = if pos_strand { 1 } else { -1 };
    let mut aa_pos = 1;
    let mut phase_index = 0;
    let mut phase: Option<i32> = None;
    let order: Vec<usize> = if pos_strand { (0..t.exons.len()).collect() } else { (0..t.exons.len()).rev().collect() };
    for i in order {
        let exon = t.exons[i];
        if let Some(prev) = prev_exon {
            let region = if pos_strand { BaseRegion::new(prev.1 + 1, exon.0 - 1) } else { BaseRegion::new(exon.1 + 1, prev.0 - 1) };
            regions.push(GeneRegionViewModel::intronic(region));
        }
        prev_exon = Some(exon);
        let (exon_start, exon_end) = (exon.0, exon.1);
        if !positions_overlap(coding_start, coding_end, exon_start, exon_end) {
            regions.push(GeneRegionViewModel::non_coding(BaseRegion::new(exon_start, exon_end)));
            continue;
        }
        let mut exon_coding_start = coding_start.max(exon_start);
        let mut exon_coding_end = coding_end.min(exon_end);
        if pos_strand {
            if exon_start < exon_coding_start {
                regions.push(GeneRegionViewModel::non_coding(BaseRegion::new(exon_start, exon_coding_start - 1)));
            }
        } else if exon_end > exon_coding_end {
            regions.push(GeneRegionViewModel::non_coding(BaseRegion::new(exon_coding_end + 1, exon_end)));
        }
        let mut pos = if pos_strand { exon_coding_start } else { exon_coding_end };
        if phase.is_none() {
            let mut p = calc_exonic_coding_phase(exon, coding_start, coding_end, pos_strand, pos);
            while p != 1 {
                p = (p + 1) % 3;
                pos += pos_inc;
                if pos_strand {
                    exon_coding_start += 1;
                } else {
                    exon_coding_end -= 1;
                }
            }
            phase = Some(p);
        }
        while position_within(pos, exon_coding_start, exon_coding_end) {
            let nuc_length = 3 - phase_index;
            let end = if pos_strand { (pos + nuc_length - 1).min(exon_coding_end) } else { (pos - nuc_length + 1).max(exon_coding_start) };
            let len = if end >= pos { end - pos + 1 } else { pos - end + 1 };
            phase_index = if len != 0 { (phase_index + len) % 3 } else { phase_index };
            regions.push(GeneRegionViewModel::amino_acid(BaseRegion::new(pos.min(end), pos.max(end)), aa_pos, b'.', b'.', false));
            if phase_index == 0 {
                aa_pos += 1;
            }
            pos = if pos_strand { end + 1 } else { end - 1 };
        }
        if pos_strand {
            if exon_end > exon_coding_end {
                regions.push(GeneRegionViewModel::non_coding(BaseRegion::new(exon_coding_end + 1, exon_end)));
            }
        } else if exon_start < exon_coding_start {
            regions.push(GeneRegionViewModel::non_coding(BaseRegion::new(exon_start, exon_coding_start - 1)));
        }
    }
    Ok(regions)
}

fn amino_acid_pos(model: &GeneRegionViewModel) -> Option<i32> {
    match model.kind {
        GeneRegionKind::AminoAcid { amino_acid_pos, .. } => Some(amino_acid_pos),
        _ => None,
    }
}

/// AminoAcidUtil.getRefGeneRegionViewModels: in genome order
pub fn ref_gene_region_view_models(t: &Transcript, gene_regions: &[GeneRegionViewModel]) -> Vec<GeneRegionViewModel> {
    let mut models = Vec::new();
    for region in gene_regions {
        match amino_acid_pos(region) {
            None => models.push(*region),
            Some(aa_pos) => {
                if aa_pos as usize > t.amino_acids.len() {
                    eprintln!(
                        "[WARN ] invalid amino-acid region({}-{}: aaIndex={}) for transcript({}:{})",
                        region.region.start, region.region.end, aa_pos, t.gene_id, t.trans_name
                    );
                    continue;
                }
                let ref_acid = t.amino_acids[aa_pos as usize - 1];
                models.push(GeneRegionViewModel::amino_acid(region.region, aa_pos, ref_acid, b'.', aa_pos == 1));
            }
        }
    }
    if !t.pos_strand {
        models.reverse();
    }
    models
}

/// Codons.codonToAminoAcid
fn codon_to_amino_acid(codon: &[u8]) -> u8 {
    const LOOKUP: &[u8; 64] = b"KNNKIIIMTTTTRSSRXYYXLFFLSSSSXCCWQHHQLLLLPPPPRRRREDDEVVVVAAAAGGGG";
    let mut index = 0;
    for b in codon {
        index = index * 4
            + match b {
                b'A' => 0,
                b'T' => 1,
                b'C' => 2,
                b'G' => 3,
                _ => return b'.',
            };
    }
    LOOKUP[index]
}

fn is_stop_codon(codon: &[u8]) -> bool {
    matches!(codon, b"TAA" | b"TAG" | b"TGA")
}

/// Nucleotides.reverseComplementBases
fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => *other,
        })
        .collect()
}

/// What the frameshift path reads from the reference genome and the variant
pub struct AltContext<'a> {
    pub variant: &'a SimpleVariant,
    pub chromosome_length: i32,
    pub ref_bases: &'a dyn Fn(i32, i32) -> Result<Vec<u8>, String>,
}

fn slice(bases: &[u8], from: i32, to: i32) -> Result<&[u8], String> {
    if from < 0 || to < from || to as usize > bases.len() {
        return Err(format!("begin {from}, end {to}, length {}", bases.len()));
    }
    Ok(&bases[from as usize..to as usize])
}

/// refBases.substring(0, variantIdx) + alt + refBases.substring(variantIdx + ref.length())
fn alt_bases_of(ref_bases: &[u8], variant: &SimpleVariant, variant_idx: i32) -> Result<Vec<u8>, String> {
    let mut alt = slice(ref_bases, 0, variant_idx)?.to_vec();
    alt.extend_from_slice(&variant.alt);
    alt.extend_from_slice(slice(ref_bases, variant_idx + variant.ref_len(), ref_bases.len() as i32)?);
    Ok(alt)
}

/// AminoAcidUtil.getAltGeneRegionViewModels: in genome order
pub fn alt_gene_region_view_models(
    t: &Transcript, gene_regions: &[GeneRegionViewModel], events: &[AminoAcidEvent], render_region: BaseRegion, alt: &AltContext,
) -> Result<Vec<GeneRegionViewModel>, String> {
    let pos_strand = t.pos_strand;
    let variant = alt.variant;
    let mut models: Vec<GeneRegionViewModel> = Vec::new();
    let mut sorted: Vec<AminoAcidEvent> = events.to_vec();
    sorted.sort_by(|x, y| (x.start_pos(), x.priority()).cmp(&(y.start_pos(), y.priority())));
    let mut event_q: VecDeque<AminoAcidEvent> = sorted.into();
    let mut stop_coding = false;
    let mut i = 0usize;
    while i < gene_regions.len() {
        let gene_region = gene_regions[i];
        if stop_coding {
            models.push(GeneRegionViewModel::non_coding(gene_region.region));
            i += 1;
            continue;
        }
        let Some(aa_pos) = amino_acid_pos(&gene_region) else {
            models.push(gene_region);
            i += 1;
            continue;
        };
        let ref_acid = *t
            .amino_acids
            .get(aa_pos as usize - 1)
            .ok_or_else(|| format!("index {} out of bounds for length {}", aa_pos - 1, t.amino_acids.len()))?;
        let mut alt_acid = ref_acid;
        while event_q.front().is_some_and(|e| e.end_pos() < aa_pos) {
            event_q.pop_front();
        }
        let mut in_del = false;
        let mut skip_rest = false;
        let mut replay: Vec<AminoAcidEvent> = Vec::new();
        while event_q.front().is_some_and(|e| e.start_pos() <= aa_pos) {
            let event = event_q.pop_front().unwrap();
            match event {
                AminoAcidEvent::Substitution { alt, .. } => {
                    replay.push(event);
                    alt_acid = alt;
                    if alt_acid == b'*' {
                        stop_coding = true;
                        break;
                    }
                }
                AminoAcidEvent::Ins { .. } => {}
                AminoAcidEvent::Del { start, end } => {
                    replay.push(event);
                    let del_len = (end - start + 1) as i64;
                    let mut j = i + 1;
                    while j < gene_regions.len() && (j as i64) < i as i64 + del_len {
                        if !gene_regions[j].is_amino_acid() {
                            break;
                        }
                        j += 1;
                    }
                    let del_end_idx = j - 1;
                    let del_start = if pos_strand { gene_region.region.start } else { gene_regions[del_end_idx].region.start };
                    let del_end = if pos_strand { gene_regions[del_end_idx].region.end } else { gene_region.region.end };
                    models.push(GeneRegionViewModel::del(BaseRegion::new(del_start, del_end)));
                    i = del_end_idx;
                    in_del = true;
                    break;
                }
                AminoAcidEvent::Ext { .. } => {
                    skip_rest = true;
                    let mut gene_region_q: VecDeque<GeneRegionViewModel> = gene_regions[i..].iter().copied().collect();
                    if pos_strand {
                        if gene_region.region.start > render_region.end {
                            break;
                        }
                        let ref_pos_start = gene_region.region.start.min(variant.pos - 1);
                        let ref_pos_end = (render_region.end + REF_BUFFER_SIZE).min(alt.chromosome_length);
                        let ref_bases = (alt.ref_bases)(ref_pos_start, ref_pos_end)?;
                        let alt_bases = alt_bases_of(&ref_bases, variant, variant.pos - ref_pos_start)?;
                        let mut idx = gene_region.region.start - ref_pos_start;
                        let mut current_aa_pos = aa_pos;
                        while idx + 2 < alt_bases.len() as i32 && idx + 2 < ref_bases.len() as i32 {
                            let region = BaseRegion::new(idx + ref_pos_start, idx + ref_pos_start + 2);
                            if region.start > render_region.end {
                                break;
                            }
                            while gene_region_q.front().is_some_and(|g| g.region.end < region.start) {
                                gene_region_q.pop_front();
                            }
                            let front = gene_region_q.front().ok_or("geneRegionQ is empty")?;
                            let mut ref_aa = b'.';
                            if front.is_amino_acid() {
                                ref_aa = codon_to_amino_acid(slice(&ref_bases, idx, idx + 3)?);
                            }
                            let codon = slice(&alt_bases, idx, idx + 3)?;
                            if is_stop_codon(codon) {
                                models.push(GeneRegionViewModel::amino_acid(region, current_aa_pos, ref_aa, b'*', false));
                                stop_coding = true;
                                break;
                            }
                            let alt_aa = codon_to_amino_acid(codon);
                            models.push(GeneRegionViewModel::amino_acid(region, current_aa_pos, ref_aa, alt_aa, current_aa_pos == 1));
                            idx += 3;
                            current_aa_pos += 1;
                        }
                        if stop_coding {
                            // S47, Java: the index into the fetched bases is used as a genome position
                            // for the non-coding remainder. Corrected: the position after the stop codon.
                            let pos = if crate::compat::java("S47") { idx + 3 } else { idx + 3 + ref_pos_start };
                            if pos <= render_region.end {
                                models.push(GeneRegionViewModel::non_coding(BaseRegion::new(pos, render_region.end)));
                            }
                        }
                    } else {
                        if gene_region.region.end < render_region.start {
                            break;
                        }
                        let ref_pos_start = (render_region.start - REF_BUFFER_SIZE).max(1);
                        let ref_pos_end = gene_region.region.end.max(variant.pos + REF_BUFFER_SIZE).min(alt.chromosome_length);
                        let ref_bases = (alt.ref_bases)(ref_pos_start, ref_pos_end)?;
                        let alt_bases = alt_bases_of(&ref_bases, variant, variant.pos - ref_pos_start)?;
                        let mut ref_idx = gene_region.region.end - ref_pos_start;
                        let mut alt_idx = ref_idx + variant.alt_len() - variant.ref_len();
                        let mut current_aa_pos = aa_pos;
                        while ref_idx - 2 >= 0 && alt_idx - 2 >= 0 {
                            let ref_region = BaseRegion::new(ref_idx + ref_pos_start - 2, ref_idx + ref_pos_start);
                            if ref_region.end < render_region.start {
                                break;
                            }
                            while gene_region_q.front().is_some_and(|g| g.region.start > ref_region.end) {
                                gene_region_q.pop_front();
                            }
                            let front = gene_region_q.front().ok_or("geneRegionQ is empty")?;
                            let mut ref_aa = b'.';
                            if front.is_amino_acid() {
                                ref_aa = codon_to_amino_acid(&reverse_complement(slice(&ref_bases, ref_idx - 2, ref_idx + 1)?));
                            }
                            let codon = reverse_complement(slice(&alt_bases, alt_idx - 2, alt_idx + 1)?);
                            if is_stop_codon(&codon) {
                                models.push(GeneRegionViewModel::amino_acid(ref_region, current_aa_pos, ref_aa, b'*', false));
                                stop_coding = true;
                                break;
                            }
                            let alt_aa = codon_to_amino_acid(&codon);
                            models.push(GeneRegionViewModel::amino_acid(ref_region, current_aa_pos, ref_aa, alt_aa, current_aa_pos == 1));
                            ref_idx -= 3;
                            alt_idx -= 3;
                            current_aa_pos += 1;
                        }
                        if stop_coding {
                            // S47, as above: an index compared with a position
                            let pos = if crate::compat::java("S47") { ref_idx - 3 } else { ref_idx - 3 + ref_pos_start };
                            if pos >= render_region.start {
                                models.push(GeneRegionViewModel::non_coding(BaseRegion::new(render_region.start, pos)));
                            }
                        }
                    }
                    break;
                }
                AminoAcidEvent::StartLost { .. } => {
                    stop_coding = true;
                    alt_acid = b'?';
                    break;
                }
            }
        }
        while let Some(event) = replay.pop() {
            event_q.push_front(event);
        }
        if skip_rest {
            break;
        }
        if in_del {
            i += 1;
            continue;
        }
        if !stop_coding || alt_acid == b'*' {
            models.push(GeneRegionViewModel::amino_acid(gene_region.region, aa_pos, ref_acid, alt_acid, aa_pos == 1));
        } else {
            models.push(GeneRegionViewModel::non_coding(gene_region.region));
        }
        i += 1;
    }
    let mut models = realign_view_models(pos_strand, variant, models);
    if !pos_strand {
        models.reverse();
    }
    Ok(models)
}

/// the region that follows `prev` in transcript direction, with the length of `model`
fn following(pos_strand: bool, prev: &GeneRegionViewModel, model: &GeneRegionViewModel) -> GeneRegionViewModel {
    let len = model.region.base_length();
    if pos_strand {
        let start = prev.region.end + 1;
        model.with_region(BaseRegion::new(start, start + len - 1))
    } else {
        let end = prev.region.start - 1;
        model.with_region(BaseRegion::new(end - len + 1, end))
    }
}

fn realign_view_models(pos_strand: bool, variant: &SimpleVariant, models: Vec<GeneRegionViewModel>) -> Vec<GeneRegionViewModel> {
    let mut realigned: Vec<GeneRegionViewModel> = Vec::new();
    for model in &models {
        if model.is_del() {
            continue;
        }
        let next = match realigned.last() {
            None => *model,
            Some(prev) => following(pos_strand, prev, model),
        };
        realigned.push(next);
    }
    let indel_size = variant.alt_len() - variant.ref_len();
    if indel_size == 0 {
        realigned
    } else if indel_size < 0 {
        realign_around_del(pos_strand, variant, realigned)
    } else {
        realign_around_ins(pos_strand, variant, realigned)
    }
}

fn realign_around_del(pos_strand: bool, variant: &SimpleVariant, models: Vec<GeneRegionViewModel>) -> Vec<GeneRegionViewModel> {
    let del_len = variant.ref_len() - variant.alt_len();
    let del_region = BaseRegion::new(variant.pos + 1, variant.pos + del_len);
    let mut realigned: Vec<GeneRegionViewModel> = Vec::new();
    let mut del_handled = false;
    let mut prev: Option<GeneRegionViewModel> = None;
    for model in &models {
        if (pos_strand && model.region.end < del_region.start) || (!pos_strand && model.region.start > del_region.end) {
            realigned.push(*model);
            prev = Some(*model);
            continue;
        }
        if del_handled {
            let next = match &prev {
                None => *model,
                Some(p) => following(pos_strand, p, model),
            };
            realigned.push(next);
            prev = realigned.last().copied();
            continue;
        }
        del_handled = true;
        if pos_strand {
            if del_region.start == model.region.start {
                let start = model.region.start;
                let end = start + del_len - 1;
                realigned.push(GeneRegionViewModel::del(BaseRegion::new(start, end)));
                let len = model.region.base_length();
                realigned.push(model.with_region(BaseRegion::new(end + 1, end + len)));
                prev = realigned.last().copied();
                continue;
            }
            let bases_before_del = del_region.start - model.region.start;
            let bases_after_del = 3 - bases_before_del;
            let start = model.region.start;
            let end = start + bases_before_del - 1;
            realigned.push(model.with_region(BaseRegion::new(start, end)));
            let start = end + 1;
            let end = start + del_len - 1;
            realigned.push(GeneRegionViewModel::del(BaseRegion::new(start, end)));
            let start = end + 1;
            let end = start + bases_after_del - 1;
            realigned.push(model.with_region(BaseRegion::new(start, end)));
        } else {
            if del_region.end == model.region.end {
                let end = model.region.end;
                let start = end - del_len + 1;
                realigned.push(GeneRegionViewModel::del(BaseRegion::new(start, end)));
                let len = model.region.base_length();
                let end = start - 1;
                realigned.push(model.with_region(BaseRegion::new(end - len + 1, end)));
                prev = realigned.last().copied();
                continue;
            }
            let bases_before_del = model.region.end - del_region.end;
            let bases_after_del = 3 - bases_before_del;
            let end = model.region.end;
            let start = end - bases_before_del + 1;
            realigned.push(model.with_region(BaseRegion::new(start, end)));
            let end = start - 1;
            let start = end - del_len + 1;
            realigned.push(GeneRegionViewModel::del(BaseRegion::new(start, end)));
            let end = start - 1;
            let start = end - bases_after_del + 1;
            realigned.push(model.with_region(BaseRegion::new(start, end)));
        }
        prev = realigned.last().copied();
    }
    realigned
}

fn realign_around_ins(pos_strand: bool, variant: &SimpleVariant, models: Vec<GeneRegionViewModel>) -> Vec<GeneRegionViewModel> {
    let ins_after_pos = variant.pos;
    let mut realigned = Vec::new();
    let mut ins_handled = false;
    for model in &models {
        if (pos_strand && model.region.end < ins_after_pos) || (!pos_strand && model.region.start > ins_after_pos + 1) || ins_handled {
            realigned.push(*model);
            continue;
        }
        ins_handled = true;
        let mut updated = *model;
        updated.bases_before_right_insert =
            if pos_strand { ins_after_pos - model.region.start + 1 } else { model.region.end - ins_after_pos };
        realigned.push(updated);
    }
    realigned
}

/// AminoAcidUtil.getGeneRegionLabel
pub fn gene_region_label(t: &Transcript, position: i32) -> String {
    let pos_strand = t.pos_strand;
    let n = t.exons.len();
    let order: Vec<usize> = if pos_strand { (0..n).collect() } else { (0..n).rev().collect() };
    let last = if pos_strand { n.wrapping_sub(1) } else { 0 };
    let mut exon_idx = 0;
    let mut in_exon = false;
    for i in order {
        let (start, end, _) = t.exons[i];
        if (pos_strand && position < start) || (!pos_strand && end < position) {
            break;
        }
        exon_idx += 1;
        if (pos_strand && position <= end) || (!pos_strand && start <= position) {
            in_exon = true;
            break;
        }
        if i == last {
            return "Downstream".to_string();
        }
    }
    if exon_idx == 0 {
        "Upstream".to_string()
    } else if in_exon {
        format!("exon {exon_idx}")
    } else {
        format!("intron {exon_idx}")
    }
}
