//! The report's chapters (orange/report/chapters), in ReportWriter's order: Front Page, Somatic Findings, Germline
//! Findings (with a reference), Immunology, RNA Findings (with ISOFOX), Tissue of Origin (with a reference),
//! Quality Control (with a QSee plot), Purity and Ploidy.

use anyhow::Result;
use hmf_plot::pdf::PageSize;

use super::common::*;
use super::tables;
use super::{DocCtx, PAGE_MARGIN_LEFT, PAGE_MARGIN_RIGHT};
use crate::compat::java;
use crate::config::OrangeConfig;
use crate::model::OrangeRecord;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Chapter {
    FrontPage,
    SomaticFindings,
    GermlineFindings,
    Immunology,
    RnaFindings,
    Cuppa,
    QualityControl,
    PurplePlots,
}

/// ReportWriter.writePdf: the chapters of a record
pub fn chapters(record: &OrangeRecord) -> Vec<Chapter> {
    let tumor_only = record.reference_id.is_none();
    let mut v = vec![Chapter::FrontPage, Chapter::SomaticFindings];
    if !tumor_only {
        v.push(Chapter::GermlineFindings);
    }
    v.push(Chapter::Immunology);
    if record.isofox.is_some() {
        v.push(Chapter::RnaFindings);
    }
    if !tumor_only {
        v.push(Chapter::Cuppa);
    }
    if record.plots.q_see_plot.is_some() {
        v.push(Chapter::QualityControl);
    }
    v.push(Chapter::PurplePlots);
    v
}

impl Chapter {
    pub fn name(&self) -> &'static str {
        match self {
            Chapter::FrontPage => "Front Page",
            Chapter::SomaticFindings => "Somatic Findings",
            Chapter::GermlineFindings => "Germline Findings",
            Chapter::Immunology => "Immunology",
            Chapter::RnaFindings => "RNA Findings",
            Chapter::Cuppa => "Tissue of Origin",
            Chapter::QualityControl => "Quality Control",
            Chapter::PurplePlots => "Purity and Ploidy",
        }
    }

    pub fn page_size(&self) -> PageSize {
        match self {
            Chapter::Cuppa | Chapter::QualityControl | Chapter::PurplePlots => PageSize::a4_landscape(),
            _ => PageSize::a4(),
        }
    }

    /// ReportChapter.contentWidth: the page width less the margins and 5
    fn content_width(&self) -> f32 {
        self.page_size().width - (5f32 + PAGE_MARGIN_LEFT + PAGE_MARGIN_RIGHT)
    }

    pub fn render(&self, ctx: &mut DocCtx, record: &OrangeRecord, config: &OrangeConfig) -> Result<()> {
        match self {
            Chapter::FrontPage => front_page(ctx, record, config),
            Chapter::SomaticFindings => somatic_findings(ctx, record, config, self.content_width()),
            Chapter::GermlineFindings => germline_findings(ctx, record, self.content_width()),
            Chapter::Immunology => immunology(ctx, record, self.content_width()),
            Chapter::RnaFindings => rna_findings(ctx, record, self.content_width()),
            Chapter::Cuppa => full_page_plot(ctx, self.name(), record.plots.cuppa_summary_plot.as_deref(), NOT_AVAILABLE),
            Chapter::QualityControl => full_page_plot(ctx, self.name(), record.plots.q_see_plot.as_deref(), "Plot not available"),
            Chapter::PurplePlots => purple_plots(ctx, record, self.content_width()),
        }
    }
}

/// addSpaceBetweenSubsections
fn space(ctx: &mut DocCtx) -> Result<()> {
    ctx.add_spacing(15.0)
}

fn has_qc_fail(record: &OrangeRecord) -> bool {
    let qc = &record.purple.fit.qc;
    tables::is_fail_no_tumor(qc) || tables::has_tumor_contaminated(qc)
}

/// FrontPageChapter.render
fn front_page(ctx: &mut DocCtx, record: &OrangeRecord, config: &OrangeConfig) -> Result<()> {
    let content_width = ctx.page_size().width - (PAGE_MARGIN_LEFT + PAGE_MARGIN_RIGHT);
    ctx.add_spacing(10.0)?;
    let t = tables::sample_summary(ctx, record, config, content_width)?;
    ctx.add_table(t)?;
    if !has_qc_fail(record) {
        ctx.add_spacing(10.0)?;
    }
    let t = tables::technical_summary(ctx, record, config, content_width)?;
    ctx.add_table(t)?;
    ctx.add_spacing(10.0)?;

    // Driver Summary and Genome Wide Biomarkers side by side
    let nested_width = (content_width / 2f32) - 3f32 * 5f32;
    let saved_y = ctx.cursor_y();
    let driver_left = PAGE_MARGIN_LEFT + 5f32 + 5f32;
    let t = tables::driver_summary(ctx, record, nested_width, driver_left)?;
    let left_final_y = ctx.add_table_no_advance(t)?;
    ctx.set_cursor_y(saved_y);
    let genome_left = ctx.margin_left() + nested_width + 4f32 * 5f32;
    let t = tables::genome_wide_features(ctx, record, nested_width, genome_left)?;
    let right_final_y = ctx.add_table_no_advance(t)?;

    let bottom_y = left_final_y.min(right_final_y) - 10f32;
    let box_width = nested_width + 10f32;
    let box_height = saved_y - bottom_y;
    ctx.draw_box_border(driver_left - 5f32, bottom_y, box_width, box_height);
    ctx.draw_box_border(genome_left - 5f32, bottom_y, box_width, box_height);
    ctx.draw_box_border(driver_left - 10f32, bottom_y - 5f32, content_width, box_height + 10f32);
    let circos_border_y = bottom_y - 5f32;
    ctx.draw_box_border(driver_left - 10f32, circos_border_y, content_width, -FRONT_CIRCOS_IMAGE_HEIGHT - 15f32);
    ctx.set_cursor_y(bottom_y - 10f32);

    let path = ctx.resolve_plot(&record.plots.purple_final_circos_plot);
    ctx.add_image(&path, content_width, FRONT_CIRCOS_IMAGE_HEIGHT)
}

/// SomaticFindingsChapter.render
fn somatic_findings(ctx: &mut DocCtx, record: &OrangeRecord, config: &OrangeConfig, width: f32) -> Result<()> {
    let res = ctx.res;
    ctx.add_paragraph("Somatic Findings", res.chapter_title_style())?;
    let tumor_only = record.reference_id.is_none();
    let purple = &record.purple;

    // small variants: the count of reported transcript impacts
    let mut count = 0;
    for v in &purple.somatic_variants {
        if v.canonical_impact.reported {
            count += 1;
        }
        count += v.other_impacts.iter().filter(|i| i.reported).count();
    }
    let t = tables::somatic_variants(ctx, &format!("Small Variants ({count})"), width, &purple.somatic_variants, tumor_only, config.rna_sample_id.is_some())?;
    ctx.add_table(t)?;
    space(ctx)?;

    let title = format!("Amplifications and Deletions ({})", purple.somatic_gains_dels.len());
    let t = tables::gains_deletions(ctx, &title, width, &purple.somatic_gains_dels, record.isofox.is_some())?;
    ctx.add_table(t)?;
    space(ctx)?;

    let title = format!("Fusions ({})", record.linx.fusions.len());
    let t = tables::dna_fusions(ctx, &title, width, &record.linx.fusions)?;
    ctx.add_table(t)?;
    space(ctx)?;

    let mut svs: Vec<i32> = record.linx.somatic_breakends.iter().map(|b| b.sv_id).collect();
    svs.sort();
    svs.dedup();
    let t = tables::disruptions(ctx, &format!("Disruptions ({})", svs.len()), width, &record.linx.somatic_breakends)?;
    ctx.add_table(t)?;
    space(ctx)?;

    if !tumor_only {
        if let Some(v) = &record.viridian {
            let title = format!("Viruses ({})", v.reportable_viruses.len());
            let t = tables::viral_presence(ctx, &title, width, &v.reportable_viruses)?;
            ctx.add_table(t)?;
        }
        space(ctx)?;
    }

    let t = tables::chr_arm_copy_numbers(ctx, "Arm Copy Number Aberrations", width, &purple.arm_copy_number_abberations)?;
    ctx.add_table(t)?;
    space(ctx)?;

    if !tumor_only {
        if let Some(sigs) = &record.sig_allocations {
            let title = format!("Signature Allocations ({})", sigs.len());
            let t = tables::signature_allocations(ctx, &title, width, sigs)?;
            ctx.add_table(t)?;
        }
        space(ctx)?;
    }

    linx_plots(ctx, record, width)
}

/// SomaticFindingsChapter.addLinxPlots: two plots per row, each with its title
fn linx_plots(ctx: &mut DocCtx, record: &OrangeRecord, width: f32) -> Result<()> {
    let res = ctx.res;
    let plots = &record.plots.linx_driver_plots;
    let col_width = width / 2f32;
    // Math.round(float) - 2
    let max_plot_height = ((width / 2f32) + 0.5f32).floor() - 2f32;
    // OR15: the plots move to a new page when the first row does not fit, the title stays; corrected: the title
    // moves with them
    let title_style = res.table_title_style();
    if !java("OR15") && !plots.is_empty() && ctx.cursor_y() - title_style.size * 1.4f32 - max_plot_height < ctx.content_end_y() {
        ctx.new_page()?;
    }
    ctx.add_paragraph(&format!("Structural Driver Plots ({})", plots.len()), title_style)?;
    let margin_left = ctx.margin_left();
    let title_style = res.table_header_style();
    let title_buffer = title_style.size * 1.4f32 * 2f32;
    let mut i = 0;
    while i < plots.len() {
        let mut row_y = ctx.cursor_y();
        if row_y - max_plot_height < ctx.content_end_y() {
            ctx.new_page()?;
            row_y = ctx.cursor_y();
        }
        let p = &plots[i];
        let path = ctx.resolve_plot(&p.filename);
        ctx.add_titled_image(&p.title, title_style, &path, margin_left, row_y, col_width, max_plot_height)?;
        if i + 1 < plots.len() {
            let p = &plots[i + 1];
            let path = ctx.resolve_plot(&p.filename);
            ctx.add_titled_image(&p.title, title_style, &path, margin_left + col_width, row_y, col_width, max_plot_height)?;
        }
        ctx.set_cursor_y(row_y - max_plot_height - 20f32 - title_buffer);
        i += 2;
    }
    Ok(())
}

/// GermlineFindingsChapter.render
fn germline_findings(ctx: &mut DocCtx, record: &OrangeRecord, width: f32) -> Result<()> {
    let res = ctx.res;
    ctx.add_paragraph("Germline Findings", res.chapter_title_style())?;
    if record.reference_id.is_none() {
        return ctx.add_paragraph(NOT_AVAILABLE, res.table_content_style());
    }
    let purple = &record.purple;
    if let (Some(_), Some(variants)) = (&purple.germline_drivers, &purple.germline_variants) {
        let t = tables::germline_variants(ctx, &format!("Small Variants ({})", variants.len()), width, variants)?;
        ctx.add_table(t)?;
    }
    space(ctx)?;
    if let Some(dels) = &purple.germline_gains_dels {
        let title = format!("Amplifications and Deletions ({})", dels.len());
        let t = tables::gains_deletions(ctx, &title, width, dels, record.isofox.is_some())?;
        ctx.add_table(t)?;
    }
    space(ctx)?;
    if let Some(b) = &record.linx.germline_breakends {
        let t = tables::disruptions(ctx, &format!("Disruptions ({})", b.len()), width, b)?;
        ctx.add_table(t)?;
    }
    space(ctx)?;
    // addGermlineCNAberrations
    let aberrations = &purple.fit.qc.germline_aberrations;
    if !aberrations.is_empty() {
        let count = aberrations.iter().filter(|a| *a != "NONE").count();
        let mut table = ctx.create_table(width)?;
        table.create_row(15.0);
        let r = table.last_row_index();
        add_cell(&res, &mut table, r, 100.0, &format!("Chromosomal Aberrations ({count})"), Style::Title)?;
        table.create_row(12.0);
        let r = table.last_row_index();
        add_cell(&res, &mut table, r, 100.0, &aberrations.join(", "), Style::Content)?;
        keep_title(&mut table, 2);
        ctx.add_table(table)?;
    }
    space(ctx)?;
    if let Some(peach) = &record.peach {
        let t = tables::pharmacogenetics(ctx, &format!("Pharmacogenetics ({})", peach.len()), width, peach)?;
        ctx.add_table(t)?;
    }
    Ok(())
}

/// ImmunologyChapter.render (the immune escape table is commented out in Java)
fn immunology(ctx: &mut DocCtx, record: &OrangeRecord, width: f32) -> Result<()> {
    let res = ctx.res;
    ctx.add_paragraph("Immunology", res.chapter_title_style())?;
    space(ctx)?;
    let Some(lilac) = &record.lilac else { return Ok(()) };
    let has_rna = record.isofox.is_some();
    let class1: Vec<_> = lilac.alleles.iter().filter(|a| a.gene_class == "MHC_CLASS_1").cloned().collect();
    let t = tables::hla_alleles(ctx, "HLA Class I Alleles", width, &class1, has_rna)?;
    ctx.add_table(t)?;
    if record.experiment_type == "WHOLE_GENOME" {
        space(ctx)?;
        let class2: Vec<_> = lilac.alleles.iter().filter(|a| a.gene_class != "MHC_CLASS_1").cloned().collect();
        let t = tables::hla_alleles(ctx, "HLA Class II Alleles", width, &class2, has_rna)?;
        ctx.add_table(t)?;
    }
    Ok(())
}

/// RnaFindingsChapter.render
fn rna_findings(ctx: &mut DocCtx, record: &OrangeRecord, width: f32) -> Result<()> {
    let res = ctx.res;
    let Some(isofox) = &record.isofox else { return Ok(()) };
    ctx.add_paragraph("RNA Findings", res.chapter_title_style())?;
    if isofox.summary.qc_status.iter().any(|s| s == "FAIL_LOW_COVERAGE") {
        ctx.add_paragraph("QC Fail", res.table_content_style())?;
    }
    let t = tables::rna_statistics(ctx, "QC", width, &isofox.summary)?;
    ctx.add_table(t)?;
    let title = format!("High Gene Expression ({})", isofox.high_expression_genes.len());
    let t = tables::expressions(ctx, &title, width, &isofox.high_expression_genes, false)?;
    ctx.add_table(t)?;
    let title = format!("Fusions Detected In RNA, Not In DNA ({})", isofox.fusions.len());
    let t = tables::rna_fusions(ctx, &title, width, &isofox.fusions)?;
    ctx.add_table(t)?;
    let title = format!("Novel Splice Junctions ({})", isofox.novel_splice_junctions.len());
    let t = tables::novel_splice_junctions(ctx, &title, width, &isofox.novel_splice_junctions)?;
    ctx.add_table(t)
}

/// CuppaChapter and QualityControlChapter: the chapter title and one plot over the page
fn full_page_plot(ctx: &mut DocCtx, name: &str, plot: Option<&str>, missing: &str) -> Result<()> {
    let res = ctx.res;
    ctx.add_paragraph(name, res.chapter_title_style())?;
    match plot {
        None => ctx.add_paragraph(missing, res.table_content_style()),
        Some(p) => {
            let path = ctx.resolve_plot(p);
            ctx.add_image(&path, FULL_PAGE_IMAGE_WIDTH, FULL_PAGE_IMAGE_HEIGHT)
        }
    }
}

/// PurplePlotsChapter.render: six plots, 3 x 2
fn purple_plots(ctx: &mut DocCtx, record: &OrangeRecord, width: f32) -> Result<()> {
    const PLOT_IMAGE_HEIGHT: f32 = 225.0;
    let res = ctx.res;
    ctx.add_paragraph("Purity and Ploidy", res.chapter_title_style())?;
    let col_width = width / 3f32;
    let margin_left = ctx.margin_left();
    let p = &record.plots;
    let mut row_y = ctx.cursor_y();
    for (k, plot) in [&p.purple_input_circos_plot, &p.purple_copy_number_plot, &p.purple_clonality_plot].iter().enumerate() {
        let path = ctx.resolve_plot(plot);
        ctx.add_image_at(&path, margin_left + k as f32 * col_width, row_y, col_width, PLOT_IMAGE_HEIGHT)?;
    }
    row_y -= PLOT_IMAGE_HEIGHT + 5f32;
    for (k, plot) in [&p.purple_purity_range_plot, &p.purple_minor_allele_map_plot, &p.purple_rainfall_plot].iter().enumerate() {
        let path = ctx.resolve_plot(plot);
        ctx.add_image_at(&path, margin_left + k as f32 * col_width, row_y, col_width, PLOT_IMAGE_HEIGHT)?;
    }
    ctx.set_cursor_y(row_y - PLOT_IMAGE_HEIGHT - 5f32);
    Ok(())
}
