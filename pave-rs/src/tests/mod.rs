//! The Java unit tests of PAVE, ported (one module per Java test class), with their test utilities
//! (ImpactTestUtils, hmf-common's MockRefGenome and GeneTestUtils.createTransExons), and tests of the port.

#![allow(unused_mut, unused_assignments, unused_variables, unused_imports, dead_code)]

/// String concatenation of the parts (Java's `a + b + c` on strings)
macro_rules! cat {
    ($($e:expr),* $(,)?) => {{
        let mut s = String::new();
        $( s.push_str(::std::convert::AsRef::<str>::as_ref(&$e)); )*
        s
    }};
}

mod del_coding;
mod hgvs_coding;
mod hgvs_coding_addresses;
mod homology_realign;
mod insert_coding;
mod mnv_coding;
mod non_coding;
mod phased_variants;
mod port;
mod protein_impact;
mod resources;
mod splice_context;
mod variant;

use std::collections::HashMap;

use anyhow::Result;

use hmf_common::ensembl::{ExonData, TranscriptData};

use crate::impact::classifier::classify_variant;
use crate::impact::{RefGenome, VariantTransImpact, NEG_STRAND, PHASE_0, PHASE_1, PHASE_2, PHASE_NONE, POS_STRAND};
use crate::variant::VariantData;

pub const CHR_1: &str = "1";
pub const GENE_ID_1: &str = "ENSG001";
pub const GENE_NAME_1: &str = "GENE_1";
pub const TRANS_ID_1: i32 = 1;
pub const CHR_2: &str = "2";
pub const CHR_3: &str = "3";
pub const GENE_ID_2: &str = "ENSG002";
pub const GENE_ID_3: &str = "ENSG003";
pub const GENE_NAME_2: &str = "GENE_2";
pub const GENE_NAME_3: &str = "GENE_3";
pub const TRANS_ID_2: i32 = 2;
pub const TRANS_ID_3: i32 = 3;

/// hmf-common's MockRefGenome: 0-based positions into each chromosome's bases, "" out of range
#[derive(Default)]
pub struct MockRefGenome {
    pub map: HashMap<String, String>,
}

impl RefGenome for MockRefGenome {
    fn base_string(&self, chromosome: &str, start: i32, end: i32) -> Result<Option<String>> {
        Ok(Some(match self.map.get(chromosome) {
            Some(b) if start >= 0 && (b.len() as i32) > end && start <= end + 1 => b[start as usize..(end + 1) as usize].to_string(),
            _ => String::new(),
        }))
    }
}

impl MockRefGenome {
    pub fn bases(&self) -> &str {
        &self.map[CHR_1]
    }
}

/// ImpactTestUtils.PAVE_DNA_BASES
pub const PAVE_DNA_BASES: [u8; 4] = [b'G', b'A', b'T', b'C'];

pub fn generate_test_bases(length: i32) -> String {
    (0..length.max(0) as usize).map(|i| PAVE_DNA_BASES[i % 4] as char).collect()
}

pub fn get_next_test_base(base: u8) -> u8 {
    match PAVE_DNA_BASES.iter().position(|b| *b == base) {
        Some(i) => PAVE_DNA_BASES[(i + 1) % 4],
        None => base,
    }
}

pub fn generate_alt(r: &str) -> String {
    r.bytes().map(|b| get_next_test_base(b) as char).collect()
}

pub fn create_mock_genome(required_bases: usize) -> MockRefGenome {
    let mut g = MockRefGenome::default();
    g.map.insert(CHR_1.to_string(), generate_test_bases(required_bases as i32));
    g
}

/// substring(start, end) of the test bases
pub fn sub(s: &str, start: i32, end: i32) -> String {
    s[start as usize..end as usize].to_string()
}

pub fn codon_bases_to_phase(mut coding_bases: i32, start_phase: i32) -> i32 {
    if start_phase == PHASE_2 {
        coding_bases += 1;
    } else if start_phase == PHASE_0 {
        coding_bases += 2;
    }
    if coding_bases >= 0 {
        return coding_bases % 3;
    }
    let m = coding_bases.abs() % 3;
    if m == 0 {
        0
    } else if m == 1 {
        2
    } else {
        1
    }
}

fn tick_phase_forward(phase: i32) -> i32 {
    (phase + 1) % 3
}

/// GeneTestUtils.createTransExons (with its negative-strand quirk of not finishing the coding region)
#[allow(clippy::too_many_arguments)]
pub fn create_trans_exons(
    gene_id: &str, trans_id: i32, strand: i8, exon_starts: &[i32], exon_length: i32, coding_start: Option<i32>, coding_end: Option<i32>,
    is_canonical: bool, biotype: &str,
) -> &'static TranscriptData {
    let exon_count = exon_starts.len();
    let trans_start = exon_starts[0];
    let trans_end = exon_starts[exon_count - 1] + exon_length;
    let mut exons: Vec<ExonData> = Vec::new();
    let has_coding_bases = coding_start.is_some() && coding_end.is_some();
    let (cs, ce) = (coding_start.unwrap_or(0), coding_end.unwrap_or(0));
    let mut in_coding = false;
    let mut finished_coding = false;
    let mut last_exon_end_phase = PHASE_NONE;

    if strand == POS_STRAND {
        for (i, &exon_start) in exon_starts.iter().enumerate() {
            let exon_end = exon_start + exon_length;
            let rank = i as i32 + 1;
            let mut exon_phase = PHASE_0;
            let mut start_phase = PHASE_NONE;
            let mut end_phase = PHASE_NONE;
            let mut exon_coding_start = 0;
            if has_coding_bases && !finished_coding {
                if !in_coding {
                    if cs <= exon_end {
                        in_coding = true;
                        exon_phase = PHASE_1;
                        start_phase = if cs == exon_start { PHASE_0 } else { PHASE_NONE };
                        exon_coding_start = cs;
                    }
                } else {
                    exon_phase = tick_phase_forward(last_exon_end_phase);
                    start_phase = last_exon_end_phase;
                    exon_coding_start = exon_start;
                }
                if in_coding {
                    let coding_bases = exon_end.min(ce) - exon_coding_start + 1;
                    end_phase = codon_bases_to_phase(coding_bases, exon_phase);
                    last_exon_end_phase = end_phase;
                    if ce <= exon_end {
                        finished_coding = true;
                        in_coding = false;
                        end_phase = PHASE_NONE;
                    }
                }
            }
            exons.push(ExonData { trans_id, start: exon_start, end: exon_end, rank, phase_start: start_phase, phase_end: end_phase });
        }
    } else {
        for i in (0..exon_count).rev() {
            let exon_start = exon_starts[i];
            let exon_end = exon_start + exon_length;
            let rank = (exon_count - i) as i32;
            let mut exon_phase = PHASE_0;
            let mut start_phase = PHASE_NONE;
            let mut end_phase = PHASE_NONE;
            let mut exon_coding_end = 0;
            if has_coding_bases && !finished_coding {
                if !in_coding {
                    if ce >= exon_start {
                        in_coding = true;
                        exon_phase = PHASE_1;
                        start_phase = if ce == exon_end { PHASE_0 } else { PHASE_NONE };
                        exon_coding_end = ce;
                    }
                } else {
                    exon_phase = tick_phase_forward(last_exon_end_phase);
                    start_phase = last_exon_end_phase;
                    exon_coding_end = exon_end;
                }
                if in_coding {
                    let coding_bases = exon_coding_end - exon_start.max(cs) + 1;
                    end_phase = codon_bases_to_phase(coding_bases, exon_phase);
                    last_exon_end_phase = end_phase;
                    if cs >= exon_start {
                        finished_coding = false;
                        in_coding = false;
                        end_phase = PHASE_NONE;
                    }
                }
            }
            exons.insert(0, ExonData { trans_id, start: exon_start, end: exon_end, rank, phase_start: start_phase, phase_end: end_phase });
        }
    }
    Box::leak(Box::new(TranscriptData {
        trans_id,
        trans_name: format!("TRAN{trans_id:04}"),
        gene_id: gene_id.to_string(),
        is_canonical,
        strand,
        trans_start,
        trans_end,
        coding_start,
        coding_end,
        bio_type: biotype.to_string(),
        ref_seq_id: None,
        exons,
    }))
}

/// codons: 15-17, 18-20, 30-32, 33-35, 36-38, 39-50, 51-53 etc
pub fn create_pos_transcript() -> &'static TranscriptData {
    create_trans_exons(GENE_ID_1, TRANS_ID_1, POS_STRAND, &[10, 30, 50, 70, 90], 10, Some(15), Some(75), false, "")
}

/// codons: 95-93, 92-90, 80-78, 77-75, 74-72, 71-60, 59-57 etc
pub fn create_neg_transcript() -> &'static TranscriptData {
    create_trans_exons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, &[10, 30, 50, 70, 90, 110, 130], 10, Some(15), Some(95), false, "")
}

pub fn create_snv(position: i32, ref_bases: &str) -> VariantData {
    let r = sub(ref_bases, position, position + 1);
    let a = (get_next_test_base(r.as_bytes()[0]) as char).to_string();
    VariantData::new(CHR_1, position, &r, &a)
}

/// AminoAcids.AMINO_ACID_TO_CODON_MAP: the codons of an amino acid
pub fn amino_acid_codons(aa: u8) -> &'static [&'static str] {
    match aa {
        b'A' => &["GCT", "GCA", "GCC", "GCG"],
        b'C' => &["TGT", "TGC"],
        b'D' => &["GAC", "GAT"],
        b'E' => &["GAG", "GAA"],
        b'F' => &["TTC", "TTT"],
        b'G' => &["GGA", "GGC", "GGG", "GGT"],
        b'H' => &["CAC", "CAT"],
        b'I' => &["ATA", "ATC", "ATT"],
        b'K' => &["AAA", "AAG"],
        b'L' => &["TTA", "TTG", "CTA", "CTC", "CTG", "CTT"],
        b'M' => &["ATG"],
        b'N' => &["AAC", "AAT"],
        b'P' => &["CCA", "CCC", "CCG", "CCT"],
        b'Q' => &["CAA", "CAG"],
        b'R' => &["AGG", "AGA", "CGA", "CGC", "CGG", "CGT"],
        b'S' => &["TCA", "AGC", "AGT", "TCC", "TCG", "TCT"],
        b'T' => &["ACA", "ACC", "ACG", "ACT"],
        b'V' => &["GTG", "GTA", "GTC", "GTT"],
        b'W' => &["TGG"],
        b'Y' => &["TAC", "TAT"],
        b'X' => &["TAA", "TAG", "TGA"],
        _ => &[],
    }
}

pub fn get_amino_acid_codon(aa: u8, index: usize) -> String {
    amino_acid_codons(aa).get(index).map(|s| s.to_string()).unwrap_or_else(|| "err".to_string())
}

pub fn get_amino_acids_codons(aas: &str, reverse_strand: bool) -> String {
    let mut out = String::new();
    for aa in aas.bytes() {
        let codon = get_amino_acid_codon(aa, 0);
        if reverse_strand {
            out = rc(&codon) + &out;
        } else {
            out += &codon;
        }
    }
    out
}

pub fn rc(s: &str) -> String {
    String::from_utf8(hmf_common::codon::reverse_complement_bases(s.as_bytes())).unwrap()
}

/// ImpactClassifier.classifyVariant for a variant in range
pub fn classify(var: &VariantData, trans: &'static TranscriptData, ref_genome: &dyn RefGenome) -> VariantTransImpact {
    classify_variant(var, trans, ref_genome).unwrap().expect("variant outside the transcript")
}

pub fn classify_opt(var: &VariantData, trans: &'static TranscriptData, ref_genome: &dyn RefGenome) -> Option<VariantTransImpact> {
    classify_variant(var, trans, ref_genome).unwrap()
}

/// AminoAcids.findAminoAcidForCodon: None for a stop or unknown codon
pub fn find_amino_acid_for_codon(codon: &str) -> Option<String> {
    if codon.len() != 3 {
        return None;
    }
    let aa = hmf_common::codon::codon_to_amino_acid(codon.as_bytes(), 0);
    if aa == hmf_common::codon::UNKNOWN || aa == b'X' {
        return None;
    }
    Some((aa as char).to_string())
}
