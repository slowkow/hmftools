//! Unit tests: NEO's Java tests (PmNeoEpitopeTest, SvNeoEpitopeTest, NeoEpitopeUtilsTest, TpmCalcsTest), ported
//! with their synthetic fixtures, and tests of the port's own (the findings, the scoring classes, the files).
//! The helpers follow hmf-common's test utilities (GeneTestUtils.createTransExons, MockRefGenome).

use anyhow::{bail, Result};
use std::collections::HashMap;

use hmf_common::ensembl::{ExonData, TranscriptData};
use hmf_common::transcript_utils::coding_bases_to_phase;

use crate::epitope::{PointMutationData, RefSource, POS_STRAND};

mod neo_epitope_utils;
mod pm_neo_epitope;
mod score;
mod sv_neo_epitope;

pub const CHR_1: &str = "1";
pub const GENE_ID_1: &str = "ENSG001";
pub const GENE_NAME_1: &str = "GENE_1";
pub const TRANS_ID_1: i32 = 1;

pub const PHASE_NONE: i32 = -1;
pub const PHASE_0: i32 = 0;
pub const PHASE_1: i32 = 1;
pub const PHASE_2: i32 = 2;

/// Java string concatenation
#[macro_export]
macro_rules! cat {
    ($($x:expr),+ $(,)?) => {{
        let mut s = String::new();
        $( s.push_str(&$x); )+
        s
    }};
}

/// String.substring(start, end)
pub fn sub(s: &str, start: i32, end: i32) -> String {
    s[start as usize..end as usize].to_string()
}

/// MockRefGenome.generateRandomBases: A, C, G, T repeated
pub fn generate_random_bases(length: usize) -> String {
    (0..length).map(|i| b"ACGT"[i % 4] as char).collect()
}

/// MockRefGenome (0-based): the bases from start to end inclusive, "" for an unknown chromosome or a range
/// outside it; Java's substring throws where the start is more than one past the end
#[derive(Default)]
pub struct MockRefGenome {
    pub bases: HashMap<String, String>,
}

impl MockRefGenome {
    pub fn put(&mut self, chromosome: &str, bases: &str) {
        self.bases.insert(chromosome.to_string(), bases.to_string());
    }
}

impl RefSource for MockRefGenome {
    fn base_string(&self, chromosome: &str, start: i32, end: i32) -> Result<String> {
        match self.bases.get(chromosome) {
            Some(b) if start >= 0 && (b.len() as i32) > end => {
                if start > end + 1 {
                    bail!("StringIndexOutOfBoundsException: begin {start}, end {}", end + 1);
                }
                Ok(b[start as usize..(end + 1) as usize].to_string())
            }
            _ => Ok(String::new()),
        }
    }
}

/// PmNeoEpitopeTest.createPointMutation
pub fn create_point_mutation(position: i32, ref_bases: &str, alt: &str) -> PointMutationData {
    PointMutationData {
        chromosome: CHR_1.to_string(),
        position,
        ref_bases: ref_bases.to_string(),
        alt: alt.to_string(),
        gene: GENE_ID_1.to_string(),
        variant_copy_number: 1.0,
        copy_number: 1.0,
        subclonal_likelihood: 0.0,
    }
}

fn tick_phase_forward(phase: i32) -> i32 {
    (phase + 1) % 3
}

/// GeneTestUtils.createTransExons (with its negative-strand quirk of not finishing the coding region); the
/// transcript is leaked so that neo-epitopes can borrow it for the test's life
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
                    end_phase = coding_bases_to_phase(coding_bases, exon_phase);
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
                    end_phase = coding_bases_to_phase(coding_bases, exon_phase);
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

/// a temporary directory for file tests, removed when dropped
pub struct TempDir(pub std::path::PathBuf);

impl TempDir {
    pub fn new(name: &str) -> TempDir {
        let p = std::env::temp_dir().join(format!("neo_rs_test_{}_{}", name, std::process::id()));
        let _ = std::fs::remove_dir_all(&p);
        std::fs::create_dir_all(&p).unwrap();
        TempDir(p)
    }
    pub fn path(&self, f: &str) -> String {
        self.0.join(f).to_string_lossy().into_owned()
    }
    pub fn dir(&self) -> String {
        format!("{}/", self.0.to_string_lossy())
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// SvNeoEpitopeTest's NeoEpitopeFusion: both breakends on CHR_1 in GENE_ID_1, copy numbers 1, chain length 0
pub fn create_fusion(
    pos_up: i32, orient_up: i8, pos_down: i32, orient_down: i8, insert_sequence: &str, transcripts: &[String; 2],
) -> hmf_common::neo::NeoEpitopeFusion {
    hmf_common::neo::NeoEpitopeFusion {
        gene_ids: [GENE_ID_1.to_string(), GENE_ID_1.to_string()],
        gene_names: [GENE_ID_1.to_string(), GENE_ID_1.to_string()],
        chromosomes: [CHR_1.to_string(), CHR_1.to_string()],
        positions: [pos_up, pos_down],
        orientations: [orient_up, orient_down],
        sv_ids: [1, 1],
        junction_copy_number: 1.0,
        copy_number: 1.0,
        insert_sequence: insert_sequence.to_string(),
        chain_length: 0,
        transcripts: transcripts.clone(),
    }
}
