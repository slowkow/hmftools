//! BreakendBuilder: forms breakends from the filtered alignments of an assembly sequence.

use crate::assembly::alignment::align_data::AlignP;
use crate::assembly::alignment::alignment_filters::filter_alignments;
use crate::assembly::alignment::alternative_alignment::AlternativeAlignment;
use crate::assembly::alignment::assembly_alignment::AsmAlignP;
use crate::assembly::alignment::breakend::{Breakend, BreakendP, BreakendSegment};
use crate::assembly::alignment::homology_data::HomologyData;
use crate::assembly::assembly_utils::reverse_complement_str;
use crate::assembly::common_utils::{positions_overlap, IndelCoords};
use crate::assembly::constants::*;
use crate::assembly::line_utils::find_line_sequence_count;
use crate::assembly::orientation::Orientation;
use crate::assembly::ref_genome::RefGenomeInterface;
use crate::assembly::sam::{get_read_index_from_position_ext, left_soft_clip_length, CigarElement, CigarOp};
use crate::assembly::types::junction_assembly::AsmP;
use crate::assembly::types::{LinkType, StructuralVariantType};

pub struct BreakendBuilder<'a> {
    ref_genome: Option<&'a dyn RefGenomeInterface>,
    assembly_alignment: AsmAlignP,
}

struct IndelBreakendData {
    chromosome: String,
    lower_position: i32,
    lower_orient: Orientation,
    upper_position: i32,
    upper_orient: Orientation,
    inserted_bases: String,
    homology: Option<HomologyData>,
    assembly_indel_start: i32,
    assembly_indel_end: i32,
}

pub fn segment_orientation(alignment: &AlignP, links_end: bool) -> Orientation {
    segment_orientation_of(alignment.orientation(), links_end)
}

pub fn segment_orientation_of(segment_orientation: Orientation, links_end: bool) -> Orientation {
    if links_end == segment_orientation.is_forward() { Orientation::Forward } else { Orientation::Reverse }
}

fn is_line_insertion(inserted_bases: &str, orientation: Orientation) -> bool {
    find_line_sequence_count(inserted_bases.as_bytes(), 0, inserted_bases.len() as i32 - 1, if orientation.is_forward() { LINE_BASE_T } else { LINE_BASE_A }) > 0
}

impl<'a> BreakendBuilder<'a> {
    pub fn new(ref_genome: Option<&'a dyn RefGenomeInterface>, assembly_alignment: AsmAlignP) -> Self {
        BreakendBuilder { ref_genome, assembly_alignment }
    }

    /// Java catches any exception here, logs it and carries on with whatever breakends were formed before it.
    /// For a SAGA-matched assembly, the breakends stated in the SAGA resource are used. False if not handled.
    pub fn form_breakends_from_saga(&self, saga_alignment: &crate::assembly::saga::SagaAlignment, align_data: AlignP) -> bool {
        // a reverse strand match should not occur and is not handled
        if !saga_alignment.is_forward() {
            return false;
        }
        let saga_assembly = &saga_alignment.saga_assembly;
        let saga_variant = &saga_assembly.variant;
        let lower_saga_breakend = &saga_variant.breakend1;
        let upper_saga_breakend = &saga_variant.breakend2;
        if lower_saga_breakend.chromosome != upper_saga_breakend.chromosome {
            return false;
        }
        let raw_insert_sequence: String;
        let indel_length: i32;
        let lower_breakend_inferred: bool;
        let upper_breakend_inferred: bool;
        if saga_variant.is_insert() {
            // the SAGA insert sequence replaces the assembled one; if the assembly did not span the whole
            // SAGA variant, the other breakend is inferred from SAGA
            raw_insert_sequence = saga_variant.insert_sequence.clone();
            indel_length = raw_insert_sequence.len() as i32;
            let saga_junction_offsets = &saga_assembly.junction_offsets;
            let start_inferred_length = saga_alignment.saga_start() - saga_junction_offsets[0];
            let end_inferred_length = saga_junction_offsets[1] - saga_alignment.saga_end();
            lower_breakend_inferred = start_inferred_length > 0;
            upper_breakend_inferred = end_inferred_length > 0;
        } else {
            raw_insert_sequence = String::new();
            indel_length = upper_saga_breakend.position - lower_saga_breakend.position - 1;
            lower_breakend_inferred = false;
            upper_breakend_inferred = false;
        }
        // the indices of the breakends within the assembly sequence, which can be out of bounds if the assembly
        // does not span the whole SAGA variant; only used for fragment support counting
        let assembly_junction_offsets = saga_alignment.query_junction_offsets();
        let assembly_indel_start = assembly_junction_offsets[0];
        let assembly_indel_end = if assembly_junction_offsets.len() > 1 { assembly_junction_offsets[1] } else { assembly_indel_start };
        let mut indel_coords = IndelCoords::new(lower_saga_breakend.position, upper_saga_breakend.position, indel_length);
        indel_coords.set_inserted_bases(raw_insert_sequence);
        let breakend_data = self.calc_indel_breakends_with_homology(&lower_saga_breakend.chromosome, &indel_coords, assembly_indel_start, assembly_indel_end);
        let (lower_breakend, upper_breakend) = self.create_indel_breakends(align_data, breakend_data);
        if lower_breakend_inferred {
            lower_breakend.m().set_saga_inferred();
        }
        if upper_breakend_inferred {
            upper_breakend.m().set_saga_inferred();
        }
        true
    }

    pub fn form_breakends(&self, alignments: &[AlignP]) {
        if let Err(e) = self.do_form_breakends(alignments) {
            crate::error!("assemblyAlign({}) error: {}", *self.assembly_alignment, e);
            for align_data in alignments {
                crate::debug!("alignment({})", **align_data);
            }
        }
    }

    fn full_sequence(&self) -> &str {
        self.assembly_alignment.full_sequence().expect("full sequence")
    }

    fn do_form_breakends(&self, alignments: &[AlignP]) -> Result<(), String> {
        if alignments.is_empty() {
            return Ok(());
        }
        let mut valid_alignments: Vec<AlignP> = Vec::new();
        let mut low_qual_alignments: Vec<AlignP> = Vec::new();
        filter_alignments(self.assembly_alignment, alignments, &mut valid_alignments, &mut low_qual_alignments);
        if valid_alignments.is_empty() {
            return Ok(());
        }
        if valid_alignments.len() == 1 {
            let single_alignment = valid_alignments[0];
            let forms_indel = self.form_indel_breakends(single_alignment);
            if !forms_indel {
                self.form_single_breakend(single_alignment, &low_qual_alignments);
            }
        } else {
            self.form_multiple_breakends(&valid_alignments, &low_qual_alignments)?;
        }
        self.set_max_local_repeat();
        Ok(())
    }

    fn is_chained(&self) -> bool {
        self.assembly_alignment.phase_set().is_some_and(|ps| ps.assembly_links().len() > 1)
    }

    fn form_indel_breakends(&self, alignment: AlignP) -> bool {
        let Some(mut indel_coords) = self.get_indel_alignment_indel_coords(alignment) else { return false };
        if indel_coords.length < MIN_INDEL_LENGTH {
            return false;
        }
        let assembly_indel_start = alignment.sequence_start() + get_read_index_from_position_ext(alignment.position_start(), alignment.cigar_elements(), indel_coords.pos_start, true, false) + 1
            - left_soft_clip_length(alignment.cigar_elements());
        let assembly_indel_end = assembly_indel_start + if indel_coords.is_insert() { indel_coords.length } else { 0 };
        if !self.is_chained() {
            let left_anchor_length = assembly_indel_start;
            let right_anchor_length = alignment.segment_length() - assembly_indel_end;
            if left_anchor_length < ALIGNMENT_INDEL_MIN_ANCHOR_LENGTH || right_anchor_length < ALIGNMENT_INDEL_MIN_ANCHOR_LENGTH {
                return false;
            }
        }
        let assembly_seq = self.full_sequence();
        if indel_coords.is_insert() && assembly_indel_start >= 0 && assembly_indel_end <= assembly_seq.len() as i32 {
            indel_coords.set_inserted_bases(assembly_seq[assembly_indel_start as usize..assembly_indel_end as usize].to_string());
        }
        let breakend_data = self.calc_indel_breakends_with_homology(alignment.chromosome(), &indel_coords, assembly_indel_start, assembly_indel_end);
        self.create_indel_breakends(alignment, breakend_data);
        true
    }

    fn get_indel_alignment_indel_coords(&self, alignment: AlignP) -> Option<IndelCoords> {
        let mut indel_coords: Option<IndelCoords> = None;
        if let Some(phase_set) = self.assembly_alignment.phase_set() {
            if phase_set.assembly_links().len() == 1 {
                let assembly_link = phase_set.assembly_links()[0];
                let mut sv_type = assembly_link.sv_type();
                if matches!(sv_type, StructuralVariantType::Del | StructuralVariantType::Dup | StructuralVariantType::Ins) {
                    let mut sv_length = assembly_link.length();
                    if sv_length == 0 && !assembly_link.inserted_bases().is_empty() {
                        sv_length = assembly_link.inserted_bases().len() as i32;
                        sv_type = StructuralVariantType::Ins;
                    }
                    let specific_indel = CigarElement::new(sv_length, if sv_type == StructuralVariantType::Del { CigarOp::D } else { CigarOp::I });
                    indel_coords = IndelCoords::find_matching_indel_coords(alignment.position_start(), alignment.cigar_elements(), &specific_indel);
                }
            }
        }
        if indel_coords.is_none() {
            indel_coords = IndelCoords::find_indel_coords(alignment.position_start(), alignment.cigar_elements(), MIN_INDEL_LENGTH);
        }
        indel_coords
    }

    fn ref_base_string(&self, chromosome: &str, start: i32, end: i32) -> String {
        self.ref_genome.expect("ref genome").get_base_string(chromosome, start, end).expect("ref genome bases")
    }

    fn calc_indel_breakends_with_homology(&self, chromosome: &str, indel_coords: &IndelCoords, mut assembly_indel_start: i32, assembly_indel_end: i32) -> IndelBreakendData {
        let mut lower_position = indel_coords.pos_start;
        let mut upper_position = indel_coords.pos_end;
        let mut lower_orient = Orientation::Forward;
        let mut upper_orient = Orientation::Reverse;
        let mut inserted_bases = indel_coords.inserted_bases().to_string();
        let mut homology: Option<HomologyData> = None;
        if indel_coords.is_insert() {
            if !inserted_bases.is_empty() {
                let hom_pos_end = upper_position;
                let bases_end = self.ref_base_string(chromosome, hom_pos_end, hom_pos_end + inserted_bases.len() as i32 - 1);
                homology = HomologyData::determine_indel_homology(&inserted_bases, &bases_end, inserted_bases.len() as i32);
                // null (unequal lengths) is a NullPointerException in Java
                let h = homology.as_ref().expect("indel homology");
                if !h.homology.is_empty() {
                    let total_inexact_homology = h.inexact_end - h.inexact_start;
                    assembly_indel_start += total_inexact_homology;
                    lower_position += 1;
                    upper_position += total_inexact_homology - 1;
                    lower_orient = Orientation::Reverse;
                    upper_orient = Orientation::Forward;
                    let mut new_homology_length = 0usize;
                    if total_inexact_homology == inserted_bases.len() as i32 {
                        let post_indel_bases = self.ref_base_string(chromosome, upper_position + 1, upper_position + inserted_bases.len() as i32);
                        let (ib, pb) = (inserted_bases.as_bytes(), post_indel_bases.as_bytes());
                        while new_homology_length < ib.len() && new_homology_length < pb.len() && ib[new_homology_length] == pb[new_homology_length] {
                            new_homology_length += 1;
                        }
                    }
                    homology = if new_homology_length > 0 {
                        Some(HomologyData::new(&inserted_bases[..new_homology_length], 0, new_homology_length as i32, 0, new_homology_length as i32))
                    } else {
                        None
                    };
                    inserted_bases = inserted_bases[total_inexact_homology as usize..].to_string();
                }
            }
        } else {
            let max_length = indel_coords.length - 1;
            let hom_pos_start = lower_position + 1;
            let bases_start = self.ref_base_string(chromosome, hom_pos_start, hom_pos_start + max_length);
            let hom_pos_end = upper_position;
            let bases_end = self.ref_base_string(chromosome, hom_pos_end, hom_pos_end + max_length);
            homology = HomologyData::determine_indel_homology(&bases_start, &bases_end, max_length);
            let h = homology.as_ref().expect("indel homology");
            if h.homology.is_empty() {
                homology = None;
            } else {
                lower_position += h.inexact_start.abs();
                upper_position += h.inexact_start.abs();
            }
        }
        IndelBreakendData { chromosome: chromosome.to_string(), lower_position, lower_orient, upper_position, upper_orient, inserted_bases, homology, assembly_indel_start, assembly_indel_end }
    }

    fn create_indel_breakends(&self, align_data: AlignP, d: IndelBreakendData) -> (BreakendP, BreakendP) {
        let aa = self.assembly_alignment;
        let lower_breakend = Breakend::new(aa, &d.chromosome, d.lower_position, d.lower_orient, &d.inserted_bases, d.homology.clone());
        let upper_breakend = Breakend::new(aa, &d.chromosome, d.upper_position, d.upper_orient, &d.inserted_bases, d.homology.clone());
        let indel_sequence_indices = Some([d.assembly_indel_start, d.assembly_indel_end]);
        lower_breakend.m().add_segment(BreakendSegment::new_indel(aa.id(), align_data.sequence_start(), 0, align_data, indel_sequence_indices));
        upper_breakend.m().add_segment(BreakendSegment::new_indel(aa.id(), align_data.sequence_start(), 0, align_data, indel_sequence_indices));
        lower_breakend.m().set_other_breakend(Some(upper_breakend));
        upper_breakend.m().set_other_breakend(Some(lower_breakend));
        aa.m().add_breakend(lower_breakend);
        aa.m().add_breakend(upper_breakend);
        (lower_breakend, upper_breakend)
    }

    fn form_single_breakend(&self, alignment: AlignP, zero_qual_alignments: &[AlignP]) {
        let full_sequence = self.full_sequence();
        let full_sequence_length = self.assembly_alignment.full_sequence_length() as usize;
        let breakend_position;
        let orientation;
        let soft_clip_length;
        if alignment.left_soft_clip_length() >= alignment.right_soft_clip_length() {
            breakend_position = alignment.position_start();
            orientation = Orientation::Reverse;
            soft_clip_length = alignment.left_soft_clip_length();
        } else {
            breakend_position = alignment.position_end();
            orientation = Orientation::Forward;
            soft_clip_length = alignment.right_soft_clip_length();
        }
        let sc = soft_clip_length as usize;
        let inserted_bases: String = if alignment.orientation().is_forward() {
            if orientation.is_reverse() { full_sequence[..sc].to_string() } else { full_sequence[full_sequence_length - sc..].to_string() }
        } else {
            let s = if orientation.is_reverse() { &full_sequence[full_sequence_length - sc..] } else { &full_sequence[..sc] };
            reverse_complement_str(s)
        };
        let required_soft_clip_length = self.calc_required_soft_clip_length(&inserted_bases, orientation);
        if soft_clip_length < required_soft_clip_length {
            return;
        }
        let aa = self.assembly_alignment;
        let breakend = Breakend::new(aa, alignment.chromosome(), breakend_position, orientation, &inserted_bases, None);
        breakend.m().add_segment(BreakendSegment::new(aa.id(), alignment.sequence_start(), 0, alignment));
        if !zero_qual_alignments.is_empty() {
            let mut alt_alignments: Vec<AlternativeAlignment> = Vec::new();
            for x in zero_qual_alignments {
                alt_alignments.extend(x.m().all_alignments());
            }
            breakend.m().set_alternative_alignments(Some(alt_alignments));
        }
        aa.m().add_breakend(breakend);
    }

    fn check_outer_single(&self, alignment: AlignP, check_start: bool, next_segment_index: i32, zero_qual_alignments: &[AlignP]) -> bool {
        let sgl_ref_base_at_end = if alignment.orientation().is_forward() { check_start } else { !check_start };
        let soft_clip_length = if sgl_ref_base_at_end { alignment.left_soft_clip_length() } else { alignment.right_soft_clip_length() };
        let full_sequence = self.full_sequence();
        let full_sequence_length = self.assembly_alignment.full_sequence_length() as usize;
        let sgl_position = if sgl_ref_base_at_end { alignment.position_start() } else { alignment.position_end() };
        let sc = soft_clip_length as usize;
        let insert_sequence = if check_start { &full_sequence[..sc] } else { &full_sequence[full_sequence_length - sc..] };
        let sgl_orientation = segment_orientation(&alignment, !check_start);
        let required_soft_clip_length = self.calc_required_soft_clip_length(insert_sequence, sgl_orientation);
        if soft_clip_length < required_soft_clip_length {
            return false;
        }
        let aa = self.assembly_alignment;
        let breakend = Breakend::new(aa, alignment.chromosome(), sgl_position, sgl_orientation, insert_sequence, None);
        breakend.m().add_segment(BreakendSegment::new(aa.id(), alignment.sequence_start(), next_segment_index, alignment));
        aa.m().add_breakend(breakend);
        let alt_alignments = if check_start { build_alt_alignments(Some(alignment), None, zero_qual_alignments) } else { build_alt_alignments(None, Some(alignment), zero_qual_alignments) };
        if !alt_alignments.is_empty() {
            breakend.m().set_alternative_alignments(Some(alt_alignments));
        }
        true
    }

    fn calc_required_soft_clip_length(&self, insert_sequence: &str, orientation: Orientation) -> i32 {
        if is_line_insertion(insert_sequence, orientation) {
            LINE_MIN_EXTENSION_LENGTH
        } else if self.assembly_alignment.is_saga_matched() {
            ALIGNMENT_MIN_SOFT_CLIP_LOWER
        } else {
            ALIGNMENT_MIN_SOFT_CLIP
        }
    }

    fn check_alignment_indel(&self, alignment: AlignP) {
        if alignment.cigar_elements().iter().any(|x| x.op.is_indel() && x.len >= MIN_INDEL_LENGTH) {
            self.form_indel_breakends(alignment);
        }
    }

    fn form_multiple_breakends(&self, alignments: &[AlignP], zero_qual_alignments: &[AlignP]) -> Result<(), String> {
        let aa = self.assembly_alignment;
        let full_sequence = self.full_sequence().to_string();
        let mut next_segment_index = 0;
        if self.check_outer_single(alignments[0], true, next_segment_index, zero_qual_alignments) {
            next_segment_index += 1;
        }
        for i in 0..alignments.len() - 1 {
            let alignment = alignments[i];
            if i == 0 {
                self.check_alignment_indel(alignment);
            }
            let (breakend_chr, mut breakend_position, breakend_orientation) = match alignment.selected_alt_alignment() {
                Some(alt) => (alt.chromosome.clone(), alt.position, alt.orient),
                None => (alignment.chromosome().to_string(), if alignment.is_forward() { alignment.position_end() } else { alignment.position_start() }, segment_orientation(&alignment, true)),
            };
            let next_alignment = alignments[i + 1];
            let (next_chr, mut next_position, next_orientation) = match next_alignment.selected_alt_alignment() {
                Some(alt) => (alt.chromosome.clone(), alt.position, alt.orient),
                None => (
                    next_alignment.chromosome().to_string(),
                    if next_alignment.is_forward() { next_alignment.position_start() } else { next_alignment.position_end() },
                    segment_orientation(&next_alignment, false),
                ),
            };
            let mut homology: Option<HomologyData> = None;
            let mut inserted_bases = String::new();
            if alignment.sequence_end() >= next_alignment.sequence_start() {
                homology = HomologyData::determine_homology(&full_sequence, &alignment, &next_alignment);
                if homology.is_none() {
                    crate::debug!("assembly({}) failed to determine homology", *aa);
                }
                if alignment.is_reverse() && next_alignment.is_reverse() {
                    match &homology {
                        Some(h) => homology = Some(h.invert(true, false)),
                        None => return Err("java.lang.NullPointerException: Cannot invoke \"com.hartwig.hmftools.esvee.assembly.alignment.HomologyData.invert(boolean, boolean)\" because \"homology\" is null".to_string()),
                    }
                }
            } else if alignment.sequence_end() < next_alignment.sequence_start() - 1 {
                let insert_seq_start = alignment.sequence_end() + 1;
                let insert_seq_end = next_alignment.sequence_start();
                inserted_bases = full_sequence[insert_seq_start as usize..insert_seq_end as usize].to_string();
            }
            let mut first_homology = homology.clone();
            let mut next_homology = homology.clone();
            if let Some(h) = &homology {
                let first_is_lower = alignment.is_lower_alignment(&next_alignment);
                let same_orientation = breakend_orientation == next_orientation;
                if !h.is_symmetrical() && same_orientation && !first_is_lower {
                    first_homology = Some(h.invert(true, false));
                }
                if same_orientation {
                    let reverse_positions = !h.is_symmetrical() && first_is_lower;
                    next_homology = Some(h.invert(reverse_positions, true));
                }
                breakend_position += first_homology.as_ref().unwrap().position_adjustment(breakend_orientation);
                next_position += next_homology.as_ref().unwrap().position_adjustment(next_orientation);
            }
            let assembly_inserted_bases = if breakend_orientation.is_forward() { inserted_bases.clone() } else { reverse_complement_str(&inserted_bases) };
            let breakend = Breakend::new(aa, &breakend_chr, breakend_position, breakend_orientation, &assembly_inserted_bases, first_homology);
            aa.m().add_breakend(breakend);
            breakend.m().add_segment(BreakendSegment::new(aa.id(), alignment.sequence_start(), next_segment_index, alignment));
            next_segment_index += 1;
            let next_inserted_bases = if next_orientation.is_reverse() { inserted_bases.clone() } else { reverse_complement_str(&inserted_bases) };
            let next_breakend = Breakend::new(aa, &next_chr, next_position, next_orientation, &next_inserted_bases, next_homology);
            aa.m().add_breakend(next_breakend);
            next_breakend.m().add_segment(BreakendSegment::new(aa.id(), next_alignment.sequence_start(), next_segment_index, next_alignment));
            next_segment_index += 1;
            breakend.m().set_other_breakend(Some(next_breakend));
            next_breakend.m().set_other_breakend(Some(breakend));
            let alt_alignments = build_alt_alignments(Some(alignment), Some(next_alignment), zero_qual_alignments);
            if !alt_alignments.is_empty() {
                breakend.m().set_alternative_alignments(Some(alt_alignments.clone()));
                next_breakend.m().set_alternative_alignments(Some(alt_alignments));
            }
            self.check_alignment_indel(next_alignment);
        }
        self.check_outer_single(alignments[alignments.len() - 1], false, next_segment_index, zero_qual_alignments);
        self.set_facing_breakends();
        Ok(())
    }

    fn set_facing_breakends(&self) {
        let aa = self.assembly_alignment;
        let breakends: Vec<BreakendP> = aa.breakends().clone();
        let n = breakends.len();
        for i in 0..n.saturating_sub(1) {
            let breakend = breakends[i];
            let next_breakend = breakends[i + 1];
            if are_facing_breakends(breakend, next_breakend) {
                breakend.m().add_facing_breakend(next_breakend);
                next_breakend.m().add_facing_breakend(breakend);
            } else if next_breakend.sv_type() == StructuralVariantType::Dup && next_breakend.orient.is_reverse() && i + 2 < n {
                let dup_other_breakend = breakends[i + 2];
                if Some(dup_other_breakend) == next_breakend.other_breakend() && are_facing_breakends(breakend, dup_other_breakend) {
                    breakend.m().add_facing_breakend(dup_other_breakend);
                    dup_other_breakend.m().add_facing_breakend(breakend);
                }
            }
        }
        let Some(phase_set) = aa.phase_set() else { return };
        if !phase_set.has_facing_links() {
            return;
        }
        for link in phase_set.assembly_links() {
            if link.type_() != LinkType::Facing {
                continue;
            }
            let first_breakend = self.find_assembly_breakend(link.first());
            let second_breakend = self.find_assembly_breakend(link.second());
            if let (Some(first_breakend), Some(second_breakend)) = (first_breakend, second_breakend) {
                if !first_breakend.facing_breakends().contains(&second_breakend) {
                    first_breakend.m().add_facing_breakend(second_breakend);
                }
                if !second_breakend.facing_breakends().contains(&first_breakend) {
                    second_breakend.m().add_facing_breakend(first_breakend);
                }
            }
        }
    }

    fn find_assembly_breakend(&self, assembly: AsmP) -> Option<BreakendP> {
        let junction = assembly.junction();
        self.assembly_alignment.breakends().iter().find(|x| x.matches_coordinates(junction.chromosome, junction.position, junction.orient)).copied()
    }

    fn set_max_local_repeat(&self) {
        let aa = self.assembly_alignment;
        if aa.breakends().len() != 2 {
            return;
        }
        if !aa.breakends().iter().all(|x| x.is_short_local_del_dup_ins()) {
            return;
        }
        let mut max_repeat_length = 0;
        for assembly in aa.assemblies() {
            let max_repeat = assembly.repeat_info().iter().map(|x| x.count).max().unwrap_or(0);
            max_repeat_length = max_repeat_length.max(max_repeat);
        }
        for x in aa.breakends() {
            x.m().set_max_local_repeat(max_repeat_length);
        }
    }
}

fn are_facing_breakends(breakend: BreakendP, next_breakend: BreakendP) -> bool {
    if breakend.other_breakend() == Some(next_breakend) {
        return false;
    }
    if breakend.orient == next_breakend.orient || breakend.chromosome != next_breakend.chromosome {
        return false;
    }
    if breakend.orient.is_reverse() {
        breakend.position < next_breakend.position && next_breakend.position - breakend.position <= PHASED_ASSEMBLY_MAX_TI
    } else {
        breakend.position > next_breakend.position && breakend.position - next_breakend.position <= PHASED_ASSEMBLY_MAX_TI
    }
}

fn build_alt_alignments(alignment_start: Option<AlignP>, alignment_end: Option<AlignP>, zero_qual_alignments: &[AlignP]) -> Vec<AlternativeAlignment> {
    let related_alignments: Vec<AlignP> = match (alignment_start, alignment_end) {
        (Some(s), Some(e)) => zero_qual_alignments.iter().filter(|x| positions_overlap(x.sequence_start(), x.sequence_end(), s.sequence_start(), e.sequence_start())).copied().collect(),
        (Some(s), None) => zero_qual_alignments.iter().filter(|x| x.sequence_start() < s.sequence_start()).copied().collect(),
        (None, Some(e)) => zero_qual_alignments.iter().filter(|x| x.sequence_start() > e.sequence_start()).copied().collect(),
        (None, None) => Vec::new(),
    };
    let mut alt_alignments: Vec<AlternativeAlignment> = Vec::new();
    for x in related_alignments {
        alt_alignments.extend(x.m().all_alignments());
    }
    alt_alignments
}
