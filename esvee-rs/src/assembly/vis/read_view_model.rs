//! ReadViewModel: one supporting read laid out against the segments of a breakend pair.

use std::collections::BTreeMap;

use crate::assembly::assembly_utils::reverse_complement_bases;
use crate::assembly::orientation::Orientation;
use crate::assembly::sam::{collapse_cigar_ops, CigarElement, CigarOp};
use crate::assembly::sequence_diff::{SequenceDiffInfo, SequenceDiffType};
use crate::assembly::types::junction_assembly::AsmP;
use crate::assembly::types::support_read::SupP;
use crate::assembly::vis::assembly_visualiser::{PairedSegmentViewModel, SegmentViewModel};
use crate::assembly::vis::{vis_cigar, VisState, READ_HEIGHT_PX};
use crate::common::vis::html::{render_read_info_table, tag, CssBuilder, CssSize, ReadInfo};
use crate::common::vis::svg::{render_base_seq, BaseRegion, BaseSeqViewModel, Canvas, BOX_PADDING};
use crate::common::vis::{VisError, VisResult};

const NO_CHROMOSOME_NAME: &str = "*";

pub struct ReadViewModel {
    support_read: SupP,
    /// index of the breakend pair in the visualiser's list
    paired_index: usize,
    read_view_model: BaseSeqViewModel,
}

fn null_pointer() -> VisError {
    VisError::new("java.lang.NullPointerException")
}

/// The operators of a cigar by the read base they belong to; operators that consume no read base are
/// listed after the base before them (key -1 at the start of the read).
fn map_cigar_operators(cigar_els: &[CigarElement]) -> BTreeMap<i32, Vec<CigarOp>> {
    let mut mapped_ops: BTreeMap<i32, Vec<CigarOp>> = BTreeMap::new();
    let mut base_idx = 0;
    for cigar_el in cigar_els {
        for _ in 0..cigar_el.len {
            if cigar_el.op.consumes_read_bases() {
                mapped_ops.entry(base_idx).or_default().push(cigar_el.op);
                base_idx += 1;
            } else {
                mapped_ops.entry(base_idx - 1).or_default().push(cigar_el.op);
            }
        }
    }
    mapped_ops
}

fn collapse_mapped_cigar_operators(mapped_ops: &BTreeMap<i32, Vec<CigarOp>>) -> Vec<CigarElement> {
    let ops: Vec<CigarOp> = mapped_ops.values().flatten().copied().collect();
    collapse_ops(&ops)
}

/// CigarUtils.collapseCigarOps, which gives no element for an empty list
fn collapse_ops(ops: &[CigarOp]) -> Vec<CigarElement> {
    if ops.is_empty() {
        return Vec::new();
    }
    collapse_cigar_ops(ops)
}

fn get_ref_view_model_index(ref_view_model: &[SegmentViewModel], read: &SupP) -> Option<usize> {
    let ref_alignment = BaseRegion::new(read.alignment_start(), read.alignment_end());
    let assembly_alignment = BaseRegion::new(read.full_assembly_index_start(), read.full_assembly_index_end());
    for (i, ref_el) in ref_view_model.iter().enumerate() {
        let Some(ref_region) = ref_el.ref_region else { continue };
        if ref_el.chromosome.as_deref() == Some(read.chromosome()) && ref_alignment.overlaps(&ref_region) && assembly_alignment.overlaps(&ref_el.assembly_region) {
            return Some(i);
        }
    }
    None
}

/// The cigar with the read's indel differences from the assembly worked in, and the shift of the read's
/// start they imply; `None` for the cigar when the read has no recorded differences (Java then keeps
/// using the read's own list).
fn correct_indels(is_built_forward: bool, mismatches: &[SequenceDiffInfo], cigar_els: &[CigarElement]) -> VisResult<(Option<Vec<CigarElement>>, i32)> {
    if mismatches.is_empty() {
        return Ok((None, 0));
    }
    let mut mapped_cigar_ops = map_cigar_operators(cigar_els);
    let mut indel_offset = 0;
    let step = if is_built_forward { 1 } else { -1 };
    for mismatch in mismatches {
        match mismatch.type_ {
            SequenceDiffType::Base => {}
            SequenceDiffType::Insert => {
                let insert_length = mismatch.indel_length.abs();
                if !is_built_forward {
                    indel_offset += insert_length;
                }
                for i in 0..insert_length {
                    if let Some(ops) = mapped_cigar_ops.get_mut(&(mismatch.read_index + i * step)) {
                        ops[0] = CigarOp::I;
                    }
                }
            }
            SequenceDiffType::Delete => {
                let del_length = mismatch.indel_length.abs();
                if !is_built_forward {
                    indel_offset -= del_length;
                }
                for _ in 0..del_length {
                    mapped_cigar_ops.get_mut(&mismatch.read_index).ok_or_else(null_pointer)?.push(CigarOp::D);
                }
            }
            SequenceDiffType::Repeat => {
                if mismatch.indel_length > 0 {
                    let insert_length = mismatch.indel_length.abs();
                    if !is_built_forward {
                        indel_offset += insert_length;
                    }
                    for i in 0..insert_length {
                        let base_idx = mismatch.read_index + if is_built_forward { 0 } else { 1 } + (i + 1) * step;
                        if let Some(ops) = mapped_cigar_ops.get_mut(&base_idx) {
                            ops[0] = CigarOp::I;
                        }
                    }
                } else {
                    let del_length = mismatch.indel_length.abs();
                    if !is_built_forward {
                        indel_offset -= del_length;
                    }
                    let base_idx = mismatch.read_index - if is_built_forward { 0 } else { 1 };
                    if let Some(ops) = mapped_cigar_ops.get_mut(&base_idx) {
                        for _ in 0..del_length {
                            ops.push(CigarOp::D);
                        }
                    }
                }
            }
            _ => return Err(VisError::new("org.apache.commons.lang3.NotImplementedException: TODO")),
        }
    }
    Ok((Some(collapse_mapped_cigar_operators(&mapped_cigar_ops)), indel_offset))
}

/// Where the two segments of the pair share bases (homology), the read is given a deletion of that
/// length after the first segment, placed where it mismatches the reference less. Java passes the read
/// bases as the base qualities here; kept.
fn handle_homology(
    ref_view_model: &PairedSegmentViewModel, read: &SupP, mut indel_offset: i32, cigar_els: &[CigarElement], read_bases: &[u8],
    read_negative_strand_flag: bool,
) -> VisResult<Option<(BaseSeqViewModel, i32)>> {
    let view_models = &ref_view_model.view_models;
    let homology_length = view_models[0].left_homology_length + view_models[view_models.len() - 1].left_homology_length;
    if homology_length <= 0 {
        return Ok(None);
    }
    let mut cigar_ops: Vec<CigarOp> = Vec::new();
    for cigar_el in cigar_els {
        for _ in 0..cigar_el.len {
            cigar_ops.push(cigar_el.op);
        }
    }
    let read_start_idx = read.full_assembly_index_start() + indel_offset + ref_view_model.read_start_offset;
    let last_idx_before_break = view_models[0].assembly_view_model.last_base_pos;
    let mut op_idx = 0;
    let mut new_cigar_ops: Vec<CigarOp> = Vec::new();
    let mut base_idx = read_start_idx;
    if base_idx <= last_idx_before_break + homology_length {
        while base_idx <= last_idx_before_break && op_idx < cigar_ops.len() {
            let current_cigar_op = cigar_ops[op_idx];
            op_idx += 1;
            new_cigar_ops.push(current_cigar_op);
            if current_cigar_op.consumes_read_bases() {
                base_idx += 1;
            }
        }
    }
    if base_idx <= last_idx_before_break {
        return Ok(None);
    }
    let read_del_length = homology_length - (base_idx - last_idx_before_break - 1);
    for _ in 0..read_del_length {
        new_cigar_ops.push(CigarOp::D);
    }
    new_cigar_ops.extend_from_slice(&cigar_ops[op_idx..]);
    let new_cigar_els = vis_cigar(&collapse_ops(&new_cigar_ops));

    let mismatch_count = |model: &BaseSeqViewModel| -> i32 { view_models.iter().filter_map(|x| x.ref_view_model.as_ref()).map(|r| model.mismatch_count(r)).sum() };
    let start = read.full_assembly_index_start() + indel_offset + ref_view_model.read_start_offset;
    let mut unshifted = BaseSeqViewModel::create(start, &new_cigar_els, read_bases, read_bases, read_negative_strand_flag)?;
    unshifted.clear_soft_clips();
    let unshifted_mismatch_count = mismatch_count(&unshifted);

    indel_offset -= read_del_length;
    let start = read.full_assembly_index_start() + indel_offset + ref_view_model.read_start_offset;
    let mut shifted = BaseSeqViewModel::create(start, &new_cigar_els, read_bases, read_bases, read_negative_strand_flag)?;
    shifted.clear_soft_clips();
    let shifted_mismatch_count = mismatch_count(&shifted);

    if unshifted_mismatch_count <= shifted_mismatch_count {
        indel_offset += read_del_length;
        return Ok(Some((unshifted, indel_offset)));
    }
    Ok(Some((shifted, indel_offset)))
}

impl ReadViewModel {
    pub fn first_base_index(&self) -> i32 {
        self.read_view_model.first_base_pos
    }
    pub fn paired_index(&self) -> usize {
        self.paired_index
    }

    /// `None` when the read lies on none of the pair's segments. As in Java, a read drawn reverse
    /// complemented has its base qualities reversed in place, and its own cigar list too when it has no
    /// recorded differences from the assembly; `state` carries both to the next use of the read.
    pub fn create(
        ref_view_model: &PairedSegmentViewModel, paired_index: usize, read: &SupP, junction_assembly: AsmP, _full_assembly_str: Option<&str>,
        state: &mut VisState,
    ) -> VisResult<Option<ReadViewModel>> {
        let mut read_negative_strand_flag = read.orientation() == Orientation::Reverse;
        let cached_read = read.cached_read().ok_or_else(null_pointer)?;
        let read_state = state.read(cached_read);

        let Some(ref_view_model_index) = get_ref_view_model_index(&ref_view_model.view_models, read) else { return Ok(None) };
        let read_ref_view_model = &ref_view_model.view_models[ref_view_model_index];
        let reverse_complemented = read_ref_view_model.is_ref_reversed != read_ref_view_model.is_assembly_reversed;
        let is_forward_junction = junction_assembly.is_forward_junction() != reverse_complemented;

        let mut indel_offset = 0;
        if is_forward_junction {
            for cigar_el in &read_state.cigar {
                if cigar_el.op == CigarOp::I {
                    indel_offset += cigar_el.len;
                } else if cigar_el.op == CigarOp::D {
                    indel_offset -= cigar_el.len;
                }
            }
        }

        let (mut corrected_cigar, index_offset) = correct_indels(ref_view_model_index == 0, cached_read.mismatches(), &read_state.cigar)?;
        indel_offset += index_offset;

        let mut read_bases: Vec<u8> = cached_read.get_bases().to_vec();
        if reverse_complemented {
            read_negative_strand_flag = !read_negative_strand_flag;
            read_bases = reverse_complement_bases(&read_bases);
            match &mut corrected_cigar {
                Some(cigar) => cigar.reverse(),
                None => read_state.cigar.reverse(),
            }
            read_state.quals.reverse();
        }
        let cigar_els: &[CigarElement] = corrected_cigar.as_deref().unwrap_or(&read_state.cigar);

        let read_view_model = match handle_homology(ref_view_model, read, indel_offset, cigar_els, &read_bases, read_negative_strand_flag)? {
            Some((model, _)) => model,
            None => {
                let start = read.full_assembly_index_start() + indel_offset + ref_view_model.read_start_offset;
                let mut model = BaseSeqViewModel::create(start, &vis_cigar(cigar_els), &read_bases, &read_state.quals, read_negative_strand_flag)?;
                model.clear_soft_clips();
                model
            }
        };
        Ok(Some(ReadViewModel { support_read: read.clone(), paired_index, read_view_model }))
    }

    /// The row of the read: its drawing and the hidden table of its properties; `None` when the read is
    /// outside the displayed part of every segment.
    pub fn render(&self, ref_view_model: &PairedSegmentViewModel) -> Option<String> {
        let read_region = BaseRegion::new(self.read_view_model.first_base_pos, self.read_view_model.last_base_pos);
        let mut total_box_width = ref_view_model.prev_box_width;
        for ref_el in &ref_view_model.view_models {
            total_box_width += ref_el.view_region.base_length() + 2 * BOX_PADDING;
        }
        let mut canvas = Canvas::new(READ_HEIGHT_PX * total_box_width as f64, READ_HEIGHT_PX);
        let init_transform = canvas.transform();
        let mut x_box_offset = ref_view_model.prev_box_width as f64;
        let mut render_count = 0;
        let count = ref_view_model.view_models.len();
        for (i, ref_el) in ref_view_model.view_models.iter().enumerate() {
            let view_region = ref_el.view_region;
            if !view_region.overlaps(&read_region) {
                x_box_offset += (view_region.base_length() + 2 * BOX_PADDING) as f64;
                continue;
            }
            let ref_bases = ref_el.ref_view_model.as_ref().unwrap_or(&ref_el.assembly_view_model);
            canvas.set_transform(init_transform);
            render_base_seq(
                &mut canvas,
                (x_box_offset, 0.0),
                READ_HEIGHT_PX,
                view_region,
                &self.read_view_model,
                true,
                Some(ref_bases),
                i == 0,
                i == count - 1,
                ref_el.is_ref_reversed,
            );
            render_count += 1;
            x_box_offset += (view_region.base_length() + 2 * BOX_PADDING) as f64;
        }
        if render_count <= 0 {
            return None;
        }

        let read = &self.support_read;
        let mate_alignment = match read.mate_chromosome() {
            Some(chromosome) if chromosome != NO_CHROMOSOME_NAME => Some((chromosome, read.mate_alignment_start(), read.mate_alignment_end())),
            _ => None,
        };
        let extra_info =
            [("Mismatch Info:", read.mismatch_info().to_string()), ("Trim count:", format!("{},{}", read.trim_count_start(), read.trim_count_end()))];
        let info = ReadInfo {
            read_name: read.id(),
            chromosome: read.chromosome(),
            alignment_start: read.alignment_start(),
            alignment_end: read.alignment_end(),
            mate_alignment,
            cigar: read.cigar(),
            insert_size: Some(read.insert_size()),
            map_q: Some(read.map_qual()),
            extra_info: &extra_info,
            ..Default::default()
        };
        let base_div_style = CssBuilder::empty().padding(CssSize::ZERO).margin(CssSize::ZERO).to_string();
        let svg_div = tag("div", &[("class", "read-svg"), ("style", &base_div_style)], &canvas.svg_element());
        Some(tag("div", &[("style", &base_div_style)], &(svg_div + &render_read_info_table(&info))))
    }
}

#[cfg(test)]
pub(crate) mod test_access {
    use super::*;

    pub fn map_and_collapse(cigar_els: &[CigarElement]) -> Vec<CigarElement> {
        collapse_mapped_cigar_operators(&map_cigar_operators(cigar_els))
    }
    pub fn correct(is_built_forward: bool, mismatches: &[SequenceDiffInfo], cigar_els: &[CigarElement]) -> VisResult<(Option<Vec<CigarElement>>, i32)> {
        correct_indels(is_built_forward, mismatches, cigar_els)
    }
}
