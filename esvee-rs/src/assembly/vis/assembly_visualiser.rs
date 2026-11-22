//! AssemblyVisualiser: the page for one assembly alignment.

use std::cmp::Ordering;
use std::path::Path;

use crate::assembly::alignment::assembly_alignment::AsmAlignP;
use crate::assembly::alignment::breakend::BreakendP;
use crate::assembly::assembly_utils::reverse_complement_str;
use crate::assembly::config::AssemblyConfig;
use crate::assembly::constants::BREAKEND_REQ_VALID_FRAGMENT_LENGTH_PERC;
use crate::assembly::orientation::Orientation;
use crate::assembly::vis::read_view_model::ReadViewModel;
use crate::assembly::vis::{vis_cigar, VisState, DISPLAY_EVERY_NTH_COORD, READ_HEIGHT_PX, VIEW_REGION_SIZE, VIS_DIR};
use crate::common::fmt::java_fixed;
use crate::common::vis::html::{base_font_style, escape, javascript, styled_table, tag, td, td_styled, CssBuilder, CssSize};
use crate::common::vis::svg::{
    render_base_seq, render_chr_labels, render_coords, string_bounds, Alignment, BaseRegion, BaseSeqViewModel, Canvas, ChrLabel, BASE_BOX_SIZE, BLACK,
    BOX_PADDING, COORD_FONT_SIZE, LIGHT_GRAY,
};
use crate::common::vis::{VisError, VisResult};

const MISSING_FIELD: &str = "missing";

/// Orientation.name()
fn orientation_name(orientation: Orientation) -> &'static str {
    if orientation.is_forward() { "FORWARD" } else { "REVERSE" }
}

fn index_error(index: i64, length: usize) -> VisError {
    VisError::new(format!("java.lang.IndexOutOfBoundsException: Index {index} out of bounds for length {length}"))
}

fn null_pointer() -> VisError {
    VisError::new("java.lang.NullPointerException")
}

/// One part of the display: the reference bases around a breakend with the assembly bases aligned to
/// them, or the inserted bases between two breakends (no reference fields).
pub struct SegmentViewModel {
    pub breakend_id: Option<i32>,
    pub chromosome: Option<String>,
    pub position: Option<i32>,
    pub ref_region: Option<BaseRegion>,
    pub view_region: BaseRegion,
    pub ref_view_region: Option<BaseRegion>,
    pub ref_view_model: Option<BaseSeqViewModel>,
    pub is_ref_reversed: bool,
    pub is_assembly_reversed: bool,
    pub assembly_view_model: BaseSeqViewModel,
    pub is_insert: bool,
    pub left_homology_length: i32,
    pub assembly_region: BaseRegion,
}

pub struct PairedSegmentViewModel {
    pub view_models: Vec<SegmentViewModel>,
    pub read_start_offset: i32,
    pub prev_box_width: i32,
}

struct PairedBreakendInfo {
    breakend: BreakendP,
    is_ref_reversed: bool,
    is_assembly_reversed: bool,
}

impl PairedBreakendInfo {
    fn new(breakend: BreakendP) -> VisResult<PairedBreakendInfo> {
        let first = breakend.segments().first().ok_or_else(|| index_error(0, 0))?;
        Ok(PairedBreakendInfo { breakend, is_ref_reversed: breakend.orient == Orientation::Reverse, is_assembly_reversed: first.alignment.orientation() == Orientation::Reverse })
    }
}

struct BreakendInfo {
    id: i32,
    chromosome: String,
    cigar_elements: Vec<crate::assembly::sam::CigarElement>,
    inserted_bases: String,
    ref_region: BaseRegion,
    assembly_seq: String,
    ref_seq: Option<String>,
    is_ref_reversed: bool,
    is_assembly_reversed: bool,
    assembly_region: BaseRegion,
    left_homology_length: i32,
}

/// java.util.PriorityQueue with a comparator: `add` each element, then `poll` until empty. Breakends that
/// compare equal come out in the order the binary heap gives, which is what decides the pairing.
fn priority_queue_order<T: Copy>(items: &[T], compare: impl Fn(&T, &T) -> Ordering) -> Vec<T> {
    let mut queue: Vec<T> = Vec::with_capacity(items.len());
    for x in items {
        // siftUpUsingComparator
        let mut k = queue.len();
        queue.push(*x);
        while k > 0 {
            let parent = (k - 1) >> 1;
            if compare(x, &queue[parent]) != Ordering::Less {
                break;
            }
            queue[k] = queue[parent];
            k = parent;
        }
        queue[k] = *x;
    }
    let mut out = Vec::with_capacity(queue.len());
    while !queue.is_empty() {
        out.push(queue[0]);
        let x = queue.pop().unwrap();
        let n = queue.len();
        if n > 0 {
            // siftDownUsingComparator
            let mut k = 0;
            let half = n >> 1;
            while k < half {
                let mut child = 2 * k + 1;
                let right = child + 1;
                if right < n && compare(&queue[child], &queue[right]) == Ordering::Greater {
                    child = right;
                }
                if compare(&x, &queue[child]) != Ordering::Greater {
                    break;
                }
                queue[k] = queue[child];
                k = child;
            }
            queue[k] = x;
        }
    }
    out
}

fn extract_breakend_info(
    config: &AssemblyConfig, state: &mut VisState, paired_breakend: &PairedBreakendInfo, is_ref_reversed: bool, full_assembly_seq: &str,
    last_seq_end: Option<i32>,
) -> VisResult<BreakendInfo> {
    let breakend = paired_breakend.breakend;
    let alignment = breakend.segments().first().ok_or_else(|| index_error(0, 0))?.alignment;
    let mut left_homology_length = 0;
    if let Some(last) = last_seq_end {
        if last >= alignment.sequence_start() {
            left_homology_length = last - alignment.sequence_start() + 1;
        }
    }
    let (begin, end) = (alignment.sequence_start(), alignment.sequence_end() + 1);
    if begin < 0 || end > full_assembly_seq.len() as i32 || begin > end {
        return Err(VisError::new(format!("java.lang.StringIndexOutOfBoundsException: begin {begin}, end {end}, length {}", full_assembly_seq.len())));
    }
    let mut assembly_seq = full_assembly_seq[begin as usize..end as usize].to_string();

    // the list of the alignment itself, or a view of it without the clip at the breakend; reversing the
    // view reverses that part of the alignment's list
    let cigar = state.align_cigar(alignment);
    let n = cigar.len();
    if n == 0 {
        return Err(index_error(if breakend.orient == Orientation::Forward { -1 } else { 0 }, 0));
    }
    let (lo, hi) = if breakend.orient == Orientation::Forward {
        if cigar[n - 1].op.is_clipping() { (0, n - 1) } else { (0, n) }
    } else if cigar[0].op.is_clipping() {
        (1, n)
    } else {
        (0, n)
    };
    let mut reverse_assembly_seq = paired_breakend.is_assembly_reversed;
    if is_ref_reversed {
        cigar[lo..hi].reverse();
        reverse_assembly_seq = !reverse_assembly_seq;
    }
    if reverse_assembly_seq {
        assembly_seq = reverse_complement_str(&assembly_seq);
        cigar[lo..hi].reverse();
    }
    let cigar_elements = cigar[lo..hi].to_vec();

    let chromosome = alignment.chromosome().to_string();
    let ref_region = BaseRegion::new(alignment.position_start(), alignment.position_end());
    let mut ref_seq = config.ref_genome().ok_or_else(null_pointer)?.get_base_string(&chromosome, ref_region.start, ref_region.end);
    if is_ref_reversed {
        ref_seq = Some(reverse_complement_str(&ref_seq.ok_or_else(null_pointer)?));
    }
    let assembly_region = BaseRegion::new(alignment.sequence_start(), alignment.sequence_end());
    Ok(BreakendInfo {
        id: breakend.id(),
        chromosome,
        cigar_elements,
        inserted_bases: breakend.inserted_bases.clone(),
        ref_region,
        assembly_seq,
        ref_seq,
        is_ref_reversed,
        is_assembly_reversed: reverse_assembly_seq,
        assembly_region,
        left_homology_length,
    })
}

/// (segments, next base index, box width)
fn get_paired_ref_view_model(
    config: &AssemblyConfig, state: &mut VisState, base_idx_start: i32, full_assembly_seq: &str, breakend_pair: &(PairedBreakendInfo, PairedBreakendInfo),
) -> VisResult<(Vec<SegmentViewModel>, i32, i32)> {
    let mut ref_view_model_pair: Vec<SegmentViewModel> = Vec::new();
    let mut last_sequence_end: Option<i32> = None;
    let mut base_idx = base_idx_start;
    let mut box_width = 0;
    for i in 0..2 {
        let paired_breakend = if i == 0 { &breakend_pair.0 } else { &breakend_pair.1 };
        let is_ref_reversed = if i == 0 { paired_breakend.is_ref_reversed } else { !paired_breakend.is_ref_reversed };
        let info = extract_breakend_info(config, state, paired_breakend, is_ref_reversed, full_assembly_seq, last_sequence_end)?;
        last_sequence_end = Some(info.assembly_region.end);

        let ref_seq_view_model = BaseSeqViewModel::from_bases(info.ref_seq.as_ref().ok_or_else(null_pointer)?.as_bytes(), base_idx);
        let assembly_seq_view_model = BaseSeqViewModel::from_string_with_cigar(info.assembly_seq.as_bytes(), &vis_cigar(&info.cigar_elements), base_idx)?;

        let view_region = if i == 0 {
            let view_region_end = base_idx + info.ref_region.base_length() - 1;
            let view_region_start = base_idx.max(view_region_end - VIEW_REGION_SIZE + 1);
            BaseRegion::new(view_region_start, view_region_end)
        } else {
            let view_region_start = base_idx;
            let view_region_end = (base_idx + info.ref_region.base_length() - 1).min(view_region_start + VIEW_REGION_SIZE - 1);
            BaseRegion::new(view_region_start, view_region_end)
        };
        base_idx += info.ref_region.base_length();

        let mut is_view_region_reversed = info.is_ref_reversed;
        if i != 0 {
            is_view_region_reversed = !is_view_region_reversed;
        }
        let (ref_view_region, pos) = if !is_view_region_reversed {
            let view_region_end = info.ref_region.end;
            let view_region_start = info.ref_region.start.max(view_region_end - VIEW_REGION_SIZE + 1);
            (BaseRegion::new(view_region_start, view_region_end), view_region_end)
        } else {
            let view_region_start = info.ref_region.start;
            let view_region_end = info.ref_region.end.min(view_region_start + VIEW_REGION_SIZE - 1);
            (BaseRegion::new(view_region_start, view_region_end), view_region_start)
        };
        ref_view_model_pair.push(SegmentViewModel {
            breakend_id: Some(info.id),
            chromosome: Some(info.chromosome.clone()),
            position: Some(pos),
            ref_region: Some(info.ref_region),
            view_region,
            ref_view_region: Some(ref_view_region),
            ref_view_model: Some(ref_seq_view_model),
            is_ref_reversed: info.is_ref_reversed,
            is_assembly_reversed: info.is_assembly_reversed,
            assembly_view_model: assembly_seq_view_model,
            is_insert: false,
            left_homology_length: info.left_homology_length,
            assembly_region: info.assembly_region,
        });
        box_width += ref_view_region.base_length() + 2 * BOX_PADDING;
        if i == 0 && !info.inserted_bases.is_empty() {
            let insert_length = info.inserted_bases.len() as i32;
            let insert_bases = if info.is_ref_reversed { reverse_complement_str(&info.inserted_bases) } else { info.inserted_bases.clone() };
            let insert_seq_view_model = BaseSeqViewModel::from_bases(insert_bases.as_bytes(), base_idx);
            let view_region = BaseRegion::new(base_idx, base_idx + insert_length - 1);
            let assembly_region = BaseRegion::new(info.assembly_region.end + 1, info.assembly_region.end + insert_length);
            base_idx += insert_length;
            ref_view_model_pair.push(SegmentViewModel {
                breakend_id: None,
                chromosome: None,
                position: None,
                ref_region: None,
                view_region,
                ref_view_region: None,
                ref_view_model: None,
                is_ref_reversed: false,
                is_assembly_reversed: false,
                assembly_view_model: insert_seq_view_model,
                is_insert: true,
                left_homology_length: 0,
                assembly_region,
            });
            box_width += insert_length + 2 * BOX_PADDING;
        }
    }
    Ok((ref_view_model_pair, base_idx, box_width))
}

fn get_ref_view_model(config: &AssemblyConfig, state: &mut VisState, assembly_alignment: AsmAlignP) -> VisResult<Vec<PairedSegmentViewModel>> {
    let breakends = assembly_alignment.breakends();
    if breakends.iter().any(|b| b.is_single()) {
        return Err(VisError::new("java.lang.RuntimeException: Got a single breakend"));
    }
    let mut ids: Vec<i32> = Vec::new();
    for breakend in breakends {
        if ids.contains(&breakend.id()) {
            crate::error!("Repeated breakend id encountered");
            std::process::exit(1);
        }
        ids.push(breakend.id());
    }
    for breakend in breakends {
        if breakend.segments().is_empty() {
            return Err(index_error(0, 0));
        }
    }
    let queue_order = priority_queue_order(breakends, |a, b| {
        let (x, y) = (a.segments()[0].alignment, b.segments()[0].alignment);
        x.sequence_start().cmp(&y.sequence_start()).then(x.sequence_end().cmp(&y.sequence_end()))
    });
    let mut breakend_pairs: Vec<(PairedBreakendInfo, PairedBreakendInfo)> = Vec::new();
    let mut seen_ids: Vec<i32> = Vec::new();
    for breakend in queue_order {
        if seen_ids.contains(&breakend.id()) {
            continue;
        }
        seen_ids.push(breakend.id());
        let mate_breakend = breakend.other_breakend().ok_or_else(null_pointer)?;
        seen_ids.push(mate_breakend.id());
        breakend_pairs.push((PairedBreakendInfo::new(breakend)?, PairedBreakendInfo::new(mate_breakend)?));
    }

    let full_assembly_seq = assembly_alignment.full_sequence();
    let mut ref_view_model: Vec<PairedSegmentViewModel> = Vec::new();
    let mut base_idx = 0;
    let mut prev_box_width = 0;
    for breakend_pair in &breakend_pairs {
        let (view_models, next_base_idx, box_width) =
            get_paired_ref_view_model(config, state, base_idx, full_assembly_seq.ok_or_else(null_pointer)?, breakend_pair)?;
        let first_breakend = breakend_pair.0.breakend;
        let read_start_offset = base_idx - first_breakend.segments()[0].alignment.sequence_start();
        ref_view_model.push(PairedSegmentViewModel { view_models, read_start_offset, prev_box_width });
        prev_box_width += box_width;
        base_idx = next_base_idx;
    }
    Ok(ref_view_model)
}

pub struct AssemblyVisualiser<'a> {
    config: &'a AssemblyConfig,
    assembly_alignment: AsmAlignP,
    ref_view_model: Vec<PairedSegmentViewModel>,
    /// reference ids, then tumor ids
    sample_names: Vec<String>,
}

type InfoRows = Vec<(String, Vec<String>)>;

/// LinkedHashMap.computeIfAbsent(key, k -> new list).add(value)
fn add_info(info: &mut InfoRows, key: &str, value: String) {
    match info.iter_mut().find(|(k, _)| k == key) {
        Some((_, values)) => values.push(value),
        None => info.push((key.to_string(), vec![value])),
    }
}

fn td_text_styled(text: &str, style: &CssBuilder) -> String {
    td_styled(&escape(text), style)
}

struct TableStyles {
    table: CssBuilder,
    header: CssBuilder,
    cell: CssBuilder,
    first_row_cell: CssBuilder,
    left_header: CssBuilder,
}

fn table_styles() -> TableStyles {
    let border = CssBuilder::empty().border(CssSize::px(1.0), "solid", BLACK).border_collapse("collapse");
    let first_row_border = CssBuilder::empty()
        .border_left(CssSize::px(1.0), "solid", BLACK)
        .border_right(CssSize::px(1.0), "solid", BLACK)
        .border_top(CssSize::px(2.0), "solid", BLACK)
        .border_bottom(CssSize::px(1.0), "solid", BLACK)
        .border_collapse("collapse");
    TableStyles {
        table: border.margin_right(CssSize::px(10.0)),
        header: CssBuilder::empty().font_weight("bold").text_align("center").background_color(LIGHT_GRAY),
        cell: CssBuilder::empty().padding_left(CssSize::em(0.5)).padding_right(CssSize::em(0.5)).merge(&border),
        first_row_cell: CssBuilder::empty().padding_left(CssSize::em(0.5)).padding_right(CssSize::em(0.5)).merge(&first_row_border),
        left_header: CssBuilder::empty()
            .background_color(LIGHT_GRAY)
            .writing_mode("vertical-rl")
            .transform("rotate(180deg)")
            .font_weight("bold")
            .text_align("center")
            .merge(&first_row_border),
    }
}

fn render_info_table_section(rows: &mut String, section_header: &str, info: &InfoRows, styles: &TableStyles) {
    let mut first_row = true;
    for (key, values) in info {
        let mut row = String::new();
        if first_row {
            row.push_str(&tag("td", &[("rowspan", &info.len().to_string()), ("style", &styles.left_header.to_string())], &escape(section_header)));
        }
        let element_style = if first_row { styles.first_row_cell.text_align("right") } else { styles.cell.text_align("right") };
        row.push_str(&td_text_styled(key, &element_style));
        for value in values {
            row.push_str(&td_text_styled(value, &element_style));
        }
        rows.push_str(&tag("tr", &[], &row));
        first_row = false;
    }
}

fn header_row(count: usize, styles: &TableStyles) -> String {
    let style = styles.cell.merge(&styles.header);
    let mut row = String::new();
    row.push_str(&td_styled("", &style));
    row.push_str(&td_styled("", &style));
    for i in 0..count {
        row.push_str(&td_text_styled(&i.to_string(), &style));
    }
    tag("tr", &[], &row)
}

impl<'a> AssemblyVisualiser<'a> {
    /// `None` when the assembly alignment has no breakends
    pub fn create(config: &'a AssemblyConfig, assembly_alignment: AsmAlignP, state: &mut VisState) -> VisResult<Option<AssemblyVisualiser<'a>>> {
        if assembly_alignment.breakends().is_empty() {
            return Ok(None);
        }
        let ref_view_model = get_ref_view_model(config, state, assembly_alignment)?;
        let sample_names: Vec<String> = config.reference_ids.iter().chain(config.tumor_ids.iter()).cloned().collect();
        Ok(Some(AssemblyVisualiser { config, assembly_alignment, ref_view_model, sample_names }))
    }

    /// index of a sample in the combined sample ids (the last one of that name, as in a HashMap filled in order)
    fn sample_name_index(&self, sample_id: &str) -> VisResult<usize> {
        self.config.combined_sample_ids().iter().rposition(|x| x == sample_id).ok_or_else(null_pointer)
    }

    pub fn write_vis_file(&mut self, state: &mut VisState) -> VisResult<()> {
        let info = self.assembly_alignment.info();
        let filename = self.get_filename(&info)?;
        let read_table_style = CssBuilder::empty().border_spacing(CssSize::ZERO);
        let vertical_spacer_style = CssBuilder::empty().height(CssSize::em(1.0));
        let vertical_spacer = tag("div", &[("style", &vertical_spacer_style.to_string())], "");
        let variant_info = self.render_variant_info()?;
        let read_table_rows = self.render_read_table(state)?;

        let top = CssBuilder::empty().vertical_align("top");
        let breakends = self.assembly_alignment.breakends().clone();
        let mut cells = String::new();
        cells.push_str(&td_styled(&self.render_variant_and_sample_support_table(&breakends)?, &top));
        cells.push_str(&td_styled(&render_properties_table(&breakends)?, &top));
        cells.push_str(&td_styled(&self.render_assembly_info(), &top));
        let sample_info = styled_table(&tag("tr", &[], &cells), &CssBuilder::empty());

        let mut body = String::with_capacity(read_table_rows.len() + (1 << 16));
        body.push_str(&variant_info);
        body.push_str(&vertical_spacer);
        body.push_str(&sample_info);
        body.push_str(&vertical_spacer);
        body.push_str(&tag("div", &[], &styled_table(&read_table_rows, &read_table_style)));
        body.push_str(&javascript());
        let mut html_str = tag("html", &[], &tag("body", &[("style", &base_font_style().to_string())], &body));
        html_str.push('\n');

        let file_path = Path::new(self.config.output_dir.as_deref().unwrap_or("")).join(VIS_DIR).join(filename);
        crate::debug!("writing assembly vis file: {}", file_path.display());
        let written = std::fs::create_dir_all(file_path.parent().unwrap()).and_then(|_| std::fs::write(&file_path, html_str));
        if let Err(e) = written {
            crate::error!("failed to write output file({}): {}", file_path.display(), e);
            std::process::exit(1);
        }
        Ok(())
    }

    fn render_variant_info(&self) -> VisResult<String> {
        let tumor = self.config.tumor_ids.first().ok_or_else(|| index_error(0, 0))?;
        let reference = self.config.reference_ids.first().ok_or_else(|| index_error(0, 0))?;
        let strong_style = CssBuilder::empty().font_weight("bold");
        Ok(tag("div", &[("style", &strong_style.to_string())], &escape(&format!("{tumor}[normal: {reference}]"))))
    }

    fn render_variant_and_sample_support_table(&self, breakends: &[BreakendP]) -> VisResult<String> {
        let styles = table_styles();
        let variant_info: InfoRows = vec![
            ("Chromosome".to_string(), breakends.iter().map(|x| x.chromosome.clone()).collect()),
            ("Position".to_string(), breakends.iter().map(|x| x.position.to_string()).collect()),
            ("Qual".to_string(), breakends.iter().map(|x| x.calc_sv_qual().to_string()).collect()),
            ("SVTYPE".to_string(), breakends.iter().map(|x| x.sv_type().name().to_string()).collect()),
        ];

        let mut sample_support_info: InfoRows = Vec::new();
        // the genotype of the first tumor sample, and the fragment totals over all samples
        let mut sample_indices: Vec<usize> = Vec::with_capacity(self.sample_names.len());
        for sample_id in &self.sample_names {
            sample_indices.push(self.sample_name_index(sample_id)?);
        }
        for breakend in breakends {
            for index in &sample_indices {
                if *index >= breakend.sample_support().len() {
                    return Err(index_error(*index as i64, breakend.sample_support().len()));
                }
            }
        }
        let tumor_sample = self.sample_names.iter().position(|x| self.config.tumor_ids.contains(x));
        for breakend in breakends {
            let mut total_split_frags = 0;
            let mut total_disc_frags = 0;
            for index in &sample_indices {
                let support = &breakend.sample_support()[*index];
                total_split_frags += support.split_fragments;
                total_disc_frags += support.discordant_fragments;
            }
            match tumor_sample {
                None => {
                    for key in ["AD", "VF", "DF", "DP", "REF", "REFPAIR", "SB", "SF"] {
                        add_info(&mut sample_support_info, key, MISSING_FIELD.to_string());
                    }
                }
                Some(sample) => {
                    let support = &breakend.sample_support()[sample_indices[sample]];
                    add_info(&mut sample_support_info, "AD", format!("0,{}", support.total_support()));
                    add_info(&mut sample_support_info, "VF", (total_disc_frags + total_split_frags).to_string());
                    add_info(&mut sample_support_info, "DF", total_disc_frags.to_string());
                    add_info(&mut sample_support_info, "DP", "0".to_string());
                    // not available in the assembly step: Java writes the literal
                    add_info(&mut sample_support_info, "REF", "TODO".to_string());
                    add_info(&mut sample_support_info, "REFPAIR", "TODO".to_string());
                    add_info(&mut sample_support_info, "SB", java_fixed(support.strand_bias(), 3));
                    add_info(&mut sample_support_info, "SF", total_split_frags.to_string());
                }
            }
        }

        let mut rows = header_row(breakends.len(), &styles);
        render_info_table_section(&mut rows, "Variant", &variant_info, &styles);
        render_info_table_section(&mut rows, "Sample support", &sample_support_info, &styles);
        Ok(tag("div", &[], &styled_table(&rows, &styles.table)))
    }

    fn render_assembly_info(&self) -> String {
        let border = CssBuilder::empty().border(CssSize::px(1.0), "solid", BLACK).border_collapse("collapse");
        let table_style = border.margin_right(CssSize::px(10.0));
        let cell_style = CssBuilder::empty().padding_left(CssSize::em(0.5)).padding_right(CssSize::em(0.5)).merge(&border).text_align("right");
        let info = [
            ("ASMINFO", self.assembly_alignment.info_limit(10)),
            ("ASMID", self.assembly_alignment.id().to_string()),
            ("ASMLEN", self.assembly_alignment.full_sequence_length().to_string()),
        ];
        let mut rows = String::new();
        for (key, value) in info {
            rows.push_str(&tag("tr", &[], &(td_text_styled(key, &cell_style) + &td_text_styled(&value, &cell_style))));
        }
        tag("div", &[], &styled_table(&rows, &table_style))
    }

    fn render_read_table(&mut self, state: &mut VisState) -> VisResult<String> {
        let light_gray_bg_style = CssBuilder::empty().background_color(LIGHT_GRAY);
        let header_style = CssBuilder::empty().background_color(LIGHT_GRAY).text_align("right");
        let mut rows = String::new();
        rows.push_str(&tag("tr", &[], &(td("") + &td_styled(&self.render_coords()?, &light_gray_bg_style))));
        let labelled = |label: &str, svg: String| tag("tr", &[], &(td_text_styled(label, &header_style) + &tag("td", &[], &svg)));
        rows.push_str(&labelled("chr", self.render_chromosome_labels()?));
        rows.push_str(&labelled("ref", self.render_ref()));
        rows.push_str(&labelled("assembly", self.render_assembly()));

        let full_assembly_str = self.assembly_alignment.full_sequence().map(|s| s.to_string());
        let mut read_view_models: Vec<ReadViewModel> = Vec::new();
        for assembly in self.assembly_alignment.assemblies() {
            for read in assembly.support().clone() {
                if read.is_supplementary() {
                    continue;
                }
                for (paired_index, paired_ref_el) in self.ref_view_model.iter().enumerate() {
                    if let Some(read_view_model) = ReadViewModel::create(paired_ref_el, paired_index, &read, assembly, full_assembly_str.as_deref(), state)? {
                        read_view_models.push(read_view_model);
                    }
                }
            }
        }
        // Collections.sort: stable
        read_view_models.sort_by_key(|x| x.first_base_index());
        for read_view_model in &read_view_models {
            let Some(read_el) = read_view_model.render(&self.ref_view_model[read_view_model.paired_index()]) else { continue };
            rows.push_str("<tr><td></td><td>");
            rows.push_str(&read_el);
            rows.push_str("</td></tr>");
        }
        Ok(rows)
    }

    fn total_box_width(&self) -> i32 {
        let mut total_box_width = 0;
        for paired_ref_el in &self.ref_view_model {
            for ref_el in &paired_ref_el.view_models {
                total_box_width += ref_el.view_region.base_length() + 2 * BOX_PADDING;
            }
        }
        total_box_width
    }

    fn render_ref(&self) -> String {
        let mut canvas = Canvas::new(READ_HEIGHT_PX * self.total_box_width() as f64, READ_HEIGHT_PX);
        let init_transform = canvas.transform();
        let mut x_box_offset = 0.0f64;
        for paired_ref_el in &self.ref_view_model {
            for ref_el in &paired_ref_el.view_models {
                let view_region = ref_el.view_region;
                if let Some(view_model) = &ref_el.ref_view_model {
                    canvas.set_transform(init_transform);
                    render_base_seq(&mut canvas, (x_box_offset, 0.0), READ_HEIGHT_PX, view_region, view_model, false, None, true, true, false);
                }
                x_box_offset += (view_region.base_length() + 2 * BOX_PADDING) as f64;
            }
        }
        canvas.svg_element()
    }

    fn render_assembly(&self) -> String {
        let mut canvas = Canvas::new(READ_HEIGHT_PX * self.total_box_width() as f64, READ_HEIGHT_PX);
        let init_transform = canvas.transform();
        let mut x_box_offset = 0.0f64;
        for paired_ref_el in &self.ref_view_model {
            for ref_el in &paired_ref_el.view_models {
                let view_region = ref_el.view_region;
                canvas.set_transform(init_transform);
                render_base_seq(
                    &mut canvas,
                    (x_box_offset, 0.0),
                    READ_HEIGHT_PX,
                    view_region,
                    &ref_el.assembly_view_model,
                    false,
                    ref_el.ref_view_model.as_ref(),
                    true,
                    true,
                    false,
                );
                x_box_offset += (view_region.base_length() + 2 * BOX_PADDING) as f64;
            }
        }
        canvas.svg_element()
    }

    fn render_coords(&self) -> VisResult<String> {
        let scaling_factor = READ_HEIGHT_PX / BASE_BOX_SIZE;
        let char_width = string_bounds(COORD_FONT_SIZE, "9").0;
        let mut max_string_width = 0.0f64;
        let mut total_box_width = 0;
        for paired_ref_el in &self.ref_view_model {
            for ref_el in &paired_ref_el.view_models {
                if ref_el.is_insert {
                    total_box_width += ref_el.view_region.base_length() + 2 * BOX_PADDING;
                    continue;
                }
                let ref_view_region = ref_el.ref_view_region.ok_or_else(null_pointer)?;
                total_box_width += ref_view_region.base_length() + 2 * BOX_PADDING;
                // the widest label: the digits all have one width, so only the length of the number matters
                for i in ref_view_region.start..=ref_view_region.end {
                    max_string_width = max_string_width.max(string_bounds(COORD_FONT_SIZE, &crate::common::vis::svg::grouped(i)).0);
                }
            }
        }
        let mut canvas =
            Canvas::new(scaling_factor * BASE_BOX_SIZE * total_box_width as f64, scaling_factor * (max_string_width + 2.0 * char_width));
        let init_transform = canvas.transform();
        let mut x_box_offset = 0.0f64;
        for paired_ref_el in &self.ref_view_model {
            for ref_el in &paired_ref_el.view_models {
                if ref_el.is_insert {
                    x_box_offset += (ref_el.view_region.base_length() + 2 * BOX_PADDING) as f64;
                    continue;
                }
                let ref_view_region = ref_el.ref_view_region.ok_or_else(null_pointer)?;
                let center_position = ref_el.position.ok_or_else(null_pointer)?;
                let canvas_size = (scaling_factor * BASE_BOX_SIZE * (ref_view_region.base_length() + 2 * BOX_PADDING) as f64, canvas.height);
                canvas.set_transform(init_transform);
                render_coords(
                    &mut canvas,
                    (x_box_offset, 0.0),
                    canvas_size,
                    READ_HEIGHT_PX,
                    ref_view_region,
                    center_position,
                    DISPLAY_EVERY_NTH_COORD,
                    max_string_width,
                    ref_el.is_ref_reversed,
                );
                x_box_offset += (ref_view_region.base_length() + 2 * BOX_PADDING) as f64;
            }
        }
        Ok(canvas.svg_element())
    }

    fn render_chromosome_labels(&self) -> VisResult<String> {
        let mut labels: Vec<ChrLabel> = Vec::new();
        let mut x_box_offset = 0;
        for paired_ref_el in &self.ref_view_model {
            for (i, ref_el) in paired_ref_el.view_models.iter().enumerate() {
                let length = ref_el.view_region.base_length() + 2 * BOX_PADDING;
                let view_region = BaseRegion::new(x_box_offset, x_box_offset + length - 1);
                if ref_el.is_insert {
                    labels.push(ChrLabel { chromosome: None, position: 0, view_region, ref_region: None, is_reversed: false, alignments: Vec::new() });
                } else {
                    labels.push(ChrLabel {
                        chromosome: ref_el.chromosome.clone(),
                        position: ref_el.position.ok_or_else(null_pointer)?,
                        view_region,
                        ref_region: ref_el.ref_region,
                        is_reversed: ref_el.is_ref_reversed,
                        alignments: vec![if i == 0 { Alignment::Right } else { Alignment::Left }],
                    });
                }
                x_box_offset += length;
            }
        }
        Ok(render_chr_labels(READ_HEIGHT_PX, &labels)?.svg_element())
    }

    fn get_filename(&self, info: &str) -> VisResult<String> {
        let sample_id = self.config.tumor_ids.first().ok_or_else(|| index_error(0, 0))?;
        Ok(format!("{}_{}.html", sample_id, info.replace(':', "_")))
    }
}

fn render_properties_table(breakends: &[BreakendP]) -> VisResult<String> {
    let styles = table_styles();
    let mut id_info: InfoRows = vec![("Id".to_string(), breakends.iter().map(|x| x.id().to_string()).collect())];
    for breakend in breakends {
        let other_breakend_id = breakend.other_breakend().ok_or_else(null_pointer)?.id();
        add_info(&mut id_info, "MATEID", other_breakend_id.to_string());
        add_info(&mut id_info, "SVID", format!("{}_{}", other_breakend_id.min(breakend.id()), other_breakend_id.max(breakend.id())));
    }

    let mut segment_properties_info: InfoRows = Vec::new();
    for breakend in breakends {
        let segments = breakend.segments();
        let join = |values: Vec<String>| values.join(",");
        add_info(&mut segment_properties_info, "BEAPOS", join(segments.iter().map(|x| x.sequence_index.to_string()).collect()));
        add_info(&mut segment_properties_info, "BEOR", orientation_name(breakend.orient).to_string());
        add_info(&mut segment_properties_info, "BEAOR", join(segments.iter().map(|x| orientation_name(x.alignment.orientation()).to_string()).collect()));
        add_info(&mut segment_properties_info, "SEGALEN", join(segments.iter().map(|x| x.alignment.aligned_bases().to_string()).collect()));
        add_info(&mut segment_properties_info, "SEGMAPQ", segments.iter().map(|x| x.alignment.map_qual()).max().unwrap_or(0).to_string());
        add_info(&mut segment_properties_info, "SEGSCO", segments.iter().map(|x| x.alignment.score()).max().unwrap_or(0).to_string());
        add_info(&mut segment_properties_info, "SEGRL", segments.iter().map(|x| x.alignment.adjusted_alignment()).max().unwrap_or(0).to_string());
    }

    let mut site_qual_info: InfoRows = Vec::new();
    for breakend in breakends {
        if breakend.is_single() && breakend.valid_fragment_length_percent() < BREAKEND_REQ_VALID_FRAGMENT_LENGTH_PERC {
            add_info(&mut site_qual_info, "AVGLEN", "0".to_string());
        } else {
            add_info(&mut site_qual_info, "AVGLEN", breakend.average_fragment_length().to_string());
        }
        add_info(&mut site_qual_info, "UFP", breakend.unique_fragment_position_count().to_string());
        add_info(&mut site_qual_info, "TPR", breakend.read_orientation_range().iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","));
        add_info(&mut site_qual_info, "MLR", breakend.max_local_repeat().to_string());
    }

    let mut rows = header_row(breakends.len(), &styles);
    render_info_table_section(&mut rows, "Id", &id_info, &styles);
    render_info_table_section(&mut rows, "Segment properties", &segment_properties_info, &styles);
    render_info_table_section(&mut rows, "Site Qual", &site_qual_info, &styles);
    Ok(tag("div", &[], &styled_table(&rows, &styles.table)))
}

#[cfg(test)]
pub(crate) fn priority_queue_order_for_test<T: Copy>(items: &[T], compare: impl Fn(&T, &T) -> Ordering) -> Vec<T> {
    priority_queue_order(items, compare)
}
