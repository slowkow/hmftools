//! BreakendBuilderSagaTest
use std::sync::Arc;

use super::*;
use crate::assembly::alignment::breakend::BreakendP;
use crate::assembly::alignment::breakend_builder::BreakendBuilder;
use crate::assembly::alignment::bwa::BwaMemAlignment;
use crate::assembly::saga::{SagaAlignment, SagaAssembly, SagaBreakend, SagaVariant};
use crate::assembly::types::StructuralVariantType::{Del, Dup, Ins};

const CHROMOSOME: &str = "chr1";
const ASSEMBLY_LENGTH: usize = 100;
const SAGA_REF_CONTEXT: usize = 100;
const SAGA_START: i32 = 50;

struct Fixture {
    _env: TestEnv,
    ref_genome: MockRefGenome,
    assembly_alignment: AsmAlignP,
}

fn fixture() -> Fixture {
    let env = test_env();
    // since the SAGA variant info is used, nothing in the assembly alignment is used in the code under test
    let assembly_alignment = AssemblyAlignment::from_assembly(create_assembly(CHROMOSOME, 1000, FORWARD, &"A".repeat(ASSEMBLY_LENGTH), ASSEMBLY_LENGTH as i32 / 2));
    Fixture { _env: env, ref_genome: MockRefGenome::new_indexing(true), assembly_alignment }
}

impl Fixture {
    fn set_up_ref_genome_for_del(&mut self, lower_breakend: i32, upper_breakend: i32, homology_length: usize) {
        let del_length = (upper_breakend - lower_breakend - 1) as usize;
        assert!(homology_length <= del_length);
        let reference = format!("{}{}{}{}", "A".repeat(lower_breakend as usize), "C".repeat(del_length), "C".repeat(homology_length), "T".repeat(del_length - homology_length));
        assert!(!self.ref_genome.ref_genome_map.contains_key(CHROMOSOME));
        self.ref_genome.put(CHROMOSOME, &reference);
    }

    fn set_up_ref_genome_for_ins(&mut self, lower_breakend: i32, insert_seq: &str, homology_length: usize) {
        assert!(homology_length <= insert_seq.len());
        // ref after the insert is homologous to the insert sequence
        let reference = format!("{}{}{}", "A".repeat(lower_breakend as usize), &insert_seq[0..homology_length], "T".repeat(insert_seq.len() - homology_length));
        assert!(!self.ref_genome.ref_genome_map.contains_key(CHROMOSOME));
        self.ref_genome.put(CHROMOSOME, &reference);
    }

    fn saga_alignment_at(&self, saga_assembly: &Arc<SagaAssembly>, saga_start: i32, cigar_str: &str) -> SagaAlignment {
        let cigar = cigar_elements_from_str(cigar_str);
        let query_length = self.assembly_alignment.full_sequence_length();
        let query_start = left_clip_length(&cigar);
        let query_end = query_length - right_clip_length(&cigar);
        let saga_end = saga_start + cigar_aligned_length(&cigar);
        let raw = BwaMemAlignment {
            sam_flag: 0,
            ref_id: 0,
            ref_start: saga_start,
            ref_end: saga_end,
            seq_start: query_start,
            seq_end: query_end,
            map_qual: 0,
            n_mismatches: 0,
            aligner_score: 0,
            suboptimal_score: 0,
            cigar: cigar_str.to_string(),
            md_tag: Some(String::new()),
            xa_tag: Some(String::new()),
            mate_ref_id: 0,
            mate_ref_start: 0,
            template_len: 0,
        };
        let saga_alignment = SagaAlignment::new(raw, cigar, query_length, saga_assembly.clone());
        saga_alignment.validate();
        saga_alignment
    }

    fn saga_alignment(&self, saga_assembly: &Arc<SagaAssembly>) -> SagaAlignment {
        self.saga_alignment_at(saga_assembly, SAGA_START, &format!("{}M", self.assembly_alignment.full_sequence_length()))
    }

    fn form_breakends(&self, saga_alignment: &SagaAlignment) -> Vec<BreakendP> {
        let align_data = P::new(AlignData::from_saga(saga_alignment));
        let builder = BreakendBuilder::new(Some(&self.ref_genome), self.assembly_alignment);
        assert!(builder.form_breakends_from_saga(saga_alignment, align_data));
        self.assembly_alignment.breakends().clone()
    }
}

fn deletion_saga_assembly(lower_pos: i32, upper_pos: i32) -> Arc<SagaAssembly> {
    let variant = SagaVariant::new("test-del", SagaBreakend::new(CHROMOSOME, lower_pos, FORWARD), SagaBreakend::new(CHROMOSOME, upper_pos, REVERSE), "").unwrap();
    let ref_start = "A".repeat(SAGA_REF_CONTEXT);
    let sequence = format!("{}{}", ref_start, "T".repeat(SAGA_REF_CONTEXT));
    Arc::new(SagaAssembly::new("", variant, vec![ref_start.len() as i32], &sequence).unwrap())
}

fn insertion_saga_assembly(lower_pos: i32, insert_sequence: &str) -> Arc<SagaAssembly> {
    let variant = SagaVariant::new("test-ins", SagaBreakend::new(CHROMOSOME, lower_pos, FORWARD), SagaBreakend::new(CHROMOSOME, lower_pos + 1, REVERSE), insert_sequence).unwrap();
    let ref_start = "A".repeat(SAGA_REF_CONTEXT);
    let ref_end = "T".repeat(SAGA_REF_CONTEXT);
    let sequence = format!("{ref_start}{insert_sequence}{ref_end}");
    let junction_offsets = vec![ref_start.len() as i32, (sequence.len() - ref_end.len()) as i32];
    Arc::new(SagaAssembly::new("", variant, junction_offsets, &sequence).unwrap())
}

fn check_breakend(breakend: BreakendP, position: i32, orientation: Orientation, sv_type: crate::assembly::types::StructuralVariantType) {
    assert_eq!(CHROMOSOME, breakend.chromosome);
    assert_eq!(position, breakend.position);
    assert_eq!(orientation, breakend.orient);
    assert_eq!(sv_type, breakend.sv_type());
}

fn indel_indices(breakend: BreakendP) -> [i32; 2] {
    breakend.segments()[0].indel_sequence_indices().expect("indel sequence indices")
}

#[test]
fn test_deletion_no_homology() {
    let mut f = fixture();
    // 200b DEL with no homology to ref
    let (lower_pos, upper_pos) = (1001, 1200);
    f.set_up_ref_genome_for_del(lower_pos, upper_pos, 0);
    let saga_assembly = deletion_saga_assembly(lower_pos, upper_pos);
    let saga_alignment = f.saga_alignment(&saga_assembly);
    let breakends = f.form_breakends(&saga_alignment);
    assert_eq!(2, breakends.len());
    check_breakend(breakends[0], lower_pos, FORWARD, Del);
    assert!(!breakends[0].is_saga_inferred());
    check_breakend(breakends[1], upper_pos, REVERSE, Del);
    assert!(!breakends[1].is_saga_inferred());
    let indices = indel_indices(breakends[0]);
    assert_eq!(saga_alignment.query_junction_offsets()[0], indices[0]);
    assert_eq!(indices[0], indices[1]);
}

#[test]
fn test_deletion_with_homology() {
    let mut f = fixture();
    // 200b DEL with 10b homology to ref
    let (lower_pos, upper_pos, homology) = (1001, 1200, 10);
    f.set_up_ref_genome_for_del(lower_pos, upper_pos, homology as usize);
    let saga_assembly = deletion_saga_assembly(lower_pos, upper_pos);
    let saga_alignment = f.saga_alignment(&saga_assembly);
    let breakends = f.form_breakends(&saga_alignment);
    assert_eq!(2, breakends.len());
    let homology_shift = homology / 2;
    check_breakend(breakends[0], lower_pos + homology_shift, FORWARD, Del);
    assert!(!breakends[0].is_saga_inferred());
    check_breakend(breakends[1], upper_pos + homology_shift, REVERSE, Del);
    assert!(!breakends[1].is_saga_inferred());
    // the indices are not shifted for deletions
    let indices = indel_indices(breakends[0]);
    assert_eq!(saga_alignment.query_junction_offsets()[0], indices[0]);
    assert_eq!(indices[0], indices[1]);
}

#[test]
fn test_insertion_no_homology() {
    let mut f = fixture();
    // 20b INS with no homology to ref
    let lower_pos = 1001;
    let insert_sequence = "ACCGTACCGTGAGAGAGAGA";
    f.set_up_ref_genome_for_ins(lower_pos, insert_sequence, 0);
    let saga_assembly = insertion_saga_assembly(lower_pos, insert_sequence);
    let saga_alignment = f.saga_alignment(&saga_assembly);
    let breakends = f.form_breakends(&saga_alignment);
    assert_eq!(2, breakends.len());
    check_breakend(breakends[0], lower_pos, FORWARD, Ins);
    assert_eq!(insert_sequence, breakends[0].inserted_bases);
    assert_eq!(CHROMOSOME, breakends[1].chromosome);
    assert_eq!(lower_pos + 1, breakends[1].position);
    assert_eq!(REVERSE, breakends[1].orient);
    let indices = indel_indices(breakends[0]);
    assert_eq!(saga_alignment.query_junction_offsets()[0], indices[0]);
    assert_eq!(saga_alignment.query_junction_offsets()[1], indices[1]);
}

#[test]
fn test_insertion_with_homology() {
    let mut f = fixture();
    // 20b INS with 10b homology to ref: converted to a DUP
    let lower_pos = 1001;
    let insert_sequence = "ACCGTTATATGAGAGAGAGA";
    let upper_pos = lower_pos + 1;
    let homology = 10;
    f.set_up_ref_genome_for_ins(lower_pos, insert_sequence, homology as usize);
    let saga_assembly = insertion_saga_assembly(lower_pos, insert_sequence);
    let saga_alignment = f.saga_alignment(&saga_assembly);
    let breakends = f.form_breakends(&saga_alignment);
    assert_eq!(2, breakends.len());
    let expected_ins_seq = &insert_sequence[homology as usize..];
    check_breakend(breakends[0], lower_pos + 1, REVERSE, Dup);
    assert_eq!(expected_ins_seq, breakends[0].inserted_bases);
    assert_eq!(CHROMOSOME, breakends[1].chromosome);
    assert_eq!(upper_pos + homology - 1, breakends[1].position);
    assert_eq!(FORWARD, breakends[1].orient);
    let indices = indel_indices(breakends[0]);
    assert_eq!(saga_alignment.query_junction_offsets()[0] + homology, indices[0]);
    assert_eq!(saga_alignment.query_junction_offsets()[1], indices[1]);
}

#[test]
fn test_insertion_inferred() {
    let mut f = fixture();
    // 500b INS which is not fully captured in the phased assembly: the second breakend is inferred
    let lower_pos = 1001;
    let insert_sequence = "C".repeat(500);
    f.set_up_ref_genome_for_ins(lower_pos, &insert_sequence, 0);
    let saga_assembly = insertion_saga_assembly(lower_pos, &insert_sequence);
    let saga_alignment = f.saga_alignment(&saga_assembly);
    let breakends = f.form_breakends(&saga_alignment);
    assert_eq!(2, breakends.len());
    check_breakend(breakends[0], lower_pos, FORWARD, Ins);
    assert_eq!(insert_sequence, breakends[0].inserted_bases);
    assert!(!breakends[0].is_saga_inferred());
    assert_eq!(CHROMOSOME, breakends[1].chromosome);
    assert_eq!(lower_pos + 1, breakends[1].position);
    assert_eq!(REVERSE, breakends[1].orient);
    assert!(breakends[1].is_saga_inferred());
    let indices = indel_indices(breakends[0]);
    assert_eq!(saga_alignment.query_junction_offsets()[0], indices[0]);
    assert_eq!(saga_alignment.query_junction_offsets()[1], indices[1]);
}

#[test]
fn test_insertion_sequence_indices_indel() {
    let mut f = fixture();
    // 30b INS with an indel when aligning to SAGA
    let lower_pos = 1001;
    let insert_sequence = "A".repeat(30);
    f.set_up_ref_genome_for_ins(lower_pos, &insert_sequence, 0);
    let saga_assembly = insertion_saga_assembly(lower_pos, &insert_sequence);
    let saga_alignment = f.saga_alignment_at(&saga_assembly, SAGA_START, "60M2I38M");
    let breakends = f.form_breakends(&saga_alignment);
    assert_eq!(2, breakends.len());
    let indices = indel_indices(breakends[0]);
    assert_eq!(saga_alignment.query_junction_offsets()[0], indices[0]);
    assert_eq!(saga_alignment.query_junction_offsets()[1], indices[1]);
}
