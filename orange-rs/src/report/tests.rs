//! Tests of the PDF report: the fonts' metrics against PDFBox's (dumped by a Java harness), the number formats
//! against Java's, boxable's line wrapping, and small reports.

use super::common::*;
use hmf_plot::pdf::{BoxMetrics, Document};

fn fonts() -> (Document, usize, usize) {
    let mut doc = Document::new();
    let r = doc.add_font(super::FONT_REGULAR.to_vec()).unwrap();
    let b = doc.add_font(super::FONT_BOLD.to_vec()).unwrap();
    (doc, r, b)
}

#[test]
fn font_metrics_equal_pdfbox() {
    // PDType0Font.load of ORANGE's fonts (FontMetricsDump): descriptor values and string widths
    let (doc, r, b) = fonts();
    let reg = doc.font(r);
    assert_eq!(reg.name(), "NimbusSanL-Regu");
    assert_eq!(reg.ascent(), 979.98047);
    assert_eq!(reg.descent(), -220.21484);
    assert_eq!(reg.x_height(), 523.9258);
    assert_eq!(reg.cap_height(), 729.0039);
    assert_eq!(reg.char_width(' ').unwrap(), 278.0);
    assert_eq!(reg.char_width('"').unwrap(), 355.0);
    assert_eq!(reg.string_width("Somatic variant:").unwrap(), 7169.0);
    assert_eq!(reg.string_width("6 (ARID1B, KRAS, NF1, SMAD4, SUZ12, TP53)").unwrap(), 21451.0);
    assert_eq!(reg.string_width("HG008-T").unwrap(), 4112.0);
    assert!((reg.average_width() - 589.0234).abs() < 0.01, "{}", reg.average_width());
    let bold = doc.font(b);
    assert_eq!(bold.name(), "NimbusSanL-Bold");
    assert_eq!(bold.x_height(), 540.03906);
    assert_eq!(bold.string_width("Somatic variant:").unwrap(), 7780.0);
    assert_eq!(bold.string_width("6 (ARID1B, KRAS, NF1, SMAD4, SUZ12, TP53)").unwrap(), 21726.0);
    assert!((bold.average_width() - 609.2529).abs() < 0.01, "{}", bold.average_width());
    // a character without a glyph fails as PDFBox does
    assert_eq!(reg.string_width("a\tb").unwrap_err().to_string(), "No glyph for U+0009 in font NimbusSanL-Regu");
    // boxable's metrics: height = (ascent - descent) / 1000 in float steps
    let m = BoxMetrics::of(reg);
    assert_eq!(m.height(7.0), (0.5239258f32 + (0.97998047f32 - 0.5239258f32) - (-0.22021484f32)) * 7.0);
}

#[test]
fn every_character_width_equals_pdfbox() {
    // the widths of the characters U+0020..U+02FF and some others, as PDFBox gives them (font_metrics.txt of the
    // harness, embedded here for the characters that occur in reports)
    let (doc, r, b) = fonts();
    let reg = [(b'A', 667.0), (b'Z', 611.0), (b'a', 556.0), (b'z', 500.0), (b'0', 556.0), (b'%', 889.0), (b'*', 389.0), (b'>', 584.0), (b'[', 278.0)];
    for (c, w) in reg {
        assert_eq!(doc.font(r).char_width(c as char).unwrap(), w, "{}", c as char);
    }
    let bold = [(b'A', 722.0), (b'a', 556.0), (b'O', 778.0), (b'R', 722.0)];
    for (c, w) in bold {
        assert_eq!(doc.font(b).char_width(c as char).unwrap(), w, "{}", c as char);
    }
}

#[test]
fn number_formats_equal_java() {
    // FormatCheck (the harness) output: value, #'%', 0'%', 0.0, 0.00, %.1f, %.0f of v*1000, %4.2e
    let cases: [(f64, &str, &str, &str, &str, &str, &str); 12] = [
        (0.0, "0%", "0.0", "0.00", "0.0", "0", "0.00e+00"),
        (-0.0, "-0%", "-0.0", "-0.00", "-0.0", "-0", "-0.00e+00"),
        (0.005, "0%", "0.0", "0.01", "0.0", "5", "5.00e-03"),
        (0.125, "12%", "0.1", "0.12", "0.1", "125", "1.25e-01"),
        (0.135, "14%", "0.1", "0.14", "0.1", "135", "1.35e-01"),
        (0.145, "14%", "0.1", "0.14", "0.1", "145", "1.45e-01"),
        (-0.004, "-0%", "-0.0", "-0.00", "-0.0", "-4", "-4.00e-03"),
        (0.9999, "100%", "1.0", "1.00", "1.0", "1000", "1.00e+00"),
        (0.57, "57%", "0.6", "0.57", "0.6", "570", "5.70e-01"),
        (1.005, "100%", "1.0", "1.00", "1.0", "1005", "1.01e+00"),
        (2.675, "268%", "2.7", "2.67", "2.7", "2675", "2.68e+00"),
        (12345.678, "1234568%", "12345.7", "12345.68", "12345.7", "12345678", "1.23e+04"),
    ];
    for (v, pct, one, two, f1, f0, e) in cases {
        assert_eq!(percentage_field(v), pct, "{v}");
        assert_eq!(percentage(v), pct, "{v}");
        assert_eq!(single_digit(v), one, "{v}");
        assert_eq!(two_digit(v), two, "{v}");
        assert_eq!(fixed(v, 1), f1, "{v}");
        assert_eq!(fixed(v * 1000.0, 0), f0, "{v}");
        assert_eq!(sci_4_2(Some(v)), e, "{v}");
    }
    assert_eq!(percentage(f64::NAN), "NaN");
    assert_eq!(sci_4_2(Some(f64::NAN)), " NaN");
    assert_eq!(sci_4_2(None), "  nu");
}

#[test]
fn zero_prefix() {
    assert_eq!(zero_prefixed("9p21.3"), "09p21.3");
    assert_eq!(zero_prefixed("chr9p21.3"), "chr9p21.3");
    assert_eq!(zero_prefixed("12q14"), "12q14");
    assert_eq!(zero_prefixed("X"), "X");
    assert_eq!(zero_prefixed("7"), "07");
}

/// The synthetic layout of the Java harness LayoutStress (work/orangepdf/harness), drawn with the port: a long scaled
/// table, a row taller than a page, side-by-side key/value tables, empty tables near the page end
fn layout_stress(ctx: &mut super::DocCtx) -> anyhow::Result<()> {
    use hmf_plot::pdf::PageSize;
    ctx.start_chapter("Layout", PageSize::a4())?;
    let width = PageSize::a4().width - 65f32;
    let res = ctx.res;
    ctx.add_paragraph("Layout", res.chapter_title_style())?;

    let headers = vec![
        HeaderData::new("Gene", 0.0, false),
        HeaderData::new("Position", 0.0, false),
        HeaderData::new("HGVS", 0.0, true),
        HeaderData::new("Value", 0.0, false),
        HeaderData::new("Flag", 0.0, false),
        HeaderData::new("Fold Change", 0.6, false),
        HeaderData::new("Driver", 0.0, false),
    ];
    let rows: Vec<Vec<String>> = (0..150).map(stress_row).collect();
    let shade: Vec<bool> = (0..150).map(|r| r % 5 == 0).collect();
    let t = scaled_table(ctx, "Long Table (150)", width, headers, rows, &shade)?;
    ctx.add_table(t)?;
    ctx.add_spacing(15.0)?;

    let mut tall = String::new();
    for i in 0..700 {
        tall.push_str(&if i % 13 == 0 { "a-much-longer-hyphenated-word".to_string() } else { format!("word{i}") });
        tall.push_str(if i % 7 == 0 { ", " } else { " " });
    }
    let rows = vec![
        vec!["tall".to_string(), tall.trim().to_string(), "end".to_string()],
        vec!["after".to_string(), "a short row after the tall one".to_string(), "x".to_string()],
    ];
    let t = fixed_table(ctx, "Tall Row", width, &[1.0, 4.0, 1.0], &["Key", "Text", "Other"], rows, &[])?;
    ctx.add_table(t)?;
    ctx.add_spacing(15.0)?;

    let nested = (width / 2f32) - 15f32;
    let saved = ctx.cursor_y();
    let mut left = ctx.create_table_at_x(nested, 40.0)?;
    key_values(&res, &mut left, "Left", 9)?;
    let left_y = ctx.add_table_no_advance(left)?;
    ctx.set_cursor_y(saved);
    let mut right = ctx.create_table_at_x(nested, 30f32 + nested + 20f32)?;
    key_values(&res, &mut right, "Right", 5)?;
    let right_y = ctx.add_table_no_advance(right)?;
    ctx.set_cursor_y(left_y.min(right_y) - 10f32);

    for i in 0..12 {
        let t = create_empty_table(ctx, &format!("Empty {i} (0)"), width)?;
        ctx.add_table(t)?;
        ctx.add_spacing(15.0)?;
    }
    Ok(())
}

fn stress_row(i: usize) -> Vec<String> {
    let mut h = String::new();
    for k in 0..=(i % 4) {
        if k > 0 {
            h.push(' ');
        }
        h.push_str(&format!("c.{}A>G [p.Lys{}Arg]", i * 3 + k, i + k));
    }
    vec![
        format!("G{i}"),
        format!("chr{}:{}", i % 22 + 1, 1000000 + i * 7919),
        h,
        if i % 37 == 0 { format!("UNBREAKABLEVALUEWITHOUTANYWRAPPOINTS{i}") } else { fixed(i as f64 / 150.0, 2) },
        (if i % 3 == 0 { "Yes" } else { "No" }).to_string(),
        if i % 11 == 0 { ">1000".to_string() } else { fixed(i as f64 * 0.37, 1) },
        (if i % 5 == 0 { "LOW" } else { "HIGH" }).to_string(),
    ]
}

fn key_values(res: &Resources, table: &mut hmf_plot::pdf::table::Table, title: &str, n: usize) -> anyhow::Result<()> {
    table.create_row(15.0);
    let r = table.last_row_index();
    add_cell(res, table, r, 100.0, title, Style::Title)?;
    for i in 0..n {
        table.create_row(COMPACT_ROW_HEIGHT);
        let r = table.last_row_index();
        add_cell(res, table, r, 50.0, &format!("Key number {i}:"), Style::Key)?;
        let v = if i % 3 == 0 { format!("a value long enough to wrap onto a second line, {i}") } else { format!("value {i}") };
        add_cell(res, table, r, 50.0, &v, Style::Value)?;
    }
    Ok(())
}

fn stress_document(java_mode: bool) -> hmf_plot::pdf::Document {
    crate::compat::with_java(java_mode, || {
        let mut ctx = super::DocCtx::new("STRESS-SAMPLE", true, None).unwrap();
        layout_stress(&mut ctx).unwrap();
        ctx.finish().unwrap()
    })
}

#[test]
fn layout_stress_pages() {
    // Java's harness (the same layout drawn with ORANGE's table classes) gives 9 pages
    let doc = stress_document(true);
    assert_eq!(doc.page_count(), 9);
    assert!(doc.to_bytes().unwrap().starts_with(b"%PDF-"));
    // written for the comparison with Java's PDF when LAYOUT_STRESS_PDF names a file
    if let Ok(path) = std::env::var("LAYOUT_STRESS_PDF") {
        doc.save(&path).unwrap();
    }
    if let Ok(path) = std::env::var("LAYOUT_STRESS_DEFAULT_PDF") {
        stress_document(false).save(&path).unwrap();
    }
}

// ------------------------------------------------------------------------------------------------
// Java's unit tests of the report classes

#[test]
fn can_resolve_plot_paths() {
    // PlotPathResolverTest
    assert_eq!(super::resolve_plot(None, "file"), "file");
    assert_eq!(super::resolve_plot(Some("/path/to/output"), "file"), "/path/to/output/file");
}

#[test]
fn scalable_table_data() {
    // TableDataTest.testScalableTableData: widths are the text lengths without ReportResources
    let len = |t: &str, _h: bool| t.chars().count() as f32;
    let (d5, d10) = ("01234".to_string(), "0123456789".to_string());
    let (d15, d20) = (format!("{d10}{d5}"), format!("{d10}{d10}"));
    let row = |a: &String, b: &String, c: &String| vec![a.clone(), b.clone(), c.clone()];
    let close = |h: &HeaderData, v: f64| assert!((h.adjusted_width() - v).abs() <= 1.0, "{} vs {v}", h.adjusted_width());

    let mut headers = vec![HeaderData::new(&d5, 0.0, false), HeaderData::new(&d5, 0.0, false), HeaderData::new(&d5, 0.0, false)];
    let rows = vec![row(&d10, &d5, &d10), row(&d5, &d10, &d10), row(&d10, &d10, &d10)];
    scale_table_widths_by(&mut headers, &rows, 30.0, len);
    close(&headers[0], 10.0);
    close(&headers[1], 10.0);
    close(&headers[2], 10.0);

    let rows = vec![row(&d15, &d5, &d10), row(&d5, &d5, &d5), row(&d10, &d5, &d5)];
    let reset = |hs: &mut Vec<HeaderData>| hs.iter_mut().for_each(|h| {
        h.header_width = 0.0;
        h.data_width = 0.0;
        h.excess_width = 0.0;
    });
    reset(&mut headers);
    scale_table_widths_by(&mut headers, &rows, 30.0, len);
    close(&headers[0], 15.0);
    close(&headers[1], 5.0);
    close(&headers[2], 10.0);

    reset(&mut headers);
    scale_table_widths_by(&mut headers, &rows, 36.0, len);
    close(&headers[0], 17.0);
    close(&headers[1], 7.0);
    close(&headers[2], 12.0);

    let mut headers = vec![HeaderData::new(&d5, 0.0, true), HeaderData::new(&d5, 0.0, false), HeaderData::new(&d5, 0.0, true)];
    let rows = vec![row(&d20, &d5, &d10), row(&d15, &d5, &d20), row(&d10, &d20, &d5)];
    scale_table_widths_by(&mut headers, &rows, 40.0, len);
    close(&headers[0], 10.0);
    close(&headers[1], 20.0);
    close(&headers[2], 10.0);

    let mut headers = vec![HeaderData::new(&d5, 0.0, false), HeaderData::new(&d20, 0.5, false), HeaderData::new(&d5, 0.0, false)];
    let rows = vec![row(&d5, &d5, &d20), row(&d5, &d5, &d5), row(&d20, &d5, &d5)];
    scale_table_widths_by(&mut headers, &rows, 50.0, len);
    close(&headers[0], 20.0);
    close(&headers[1], 10.0);
    close(&headers[2], 20.0);
}

#[test]
fn can_sort_signature_allocations() {
    // SignatureAllocationTableTest
    let create = |s: &str, a: f64| crate::model::SignatureAllocation { signature: s.into(), etiology: "test".into(), allocation: a, percent: 0.0 };
    let (a1, a2, a3, a4) = (create("Sig1", 120.0), create("Sig2", 200.0), create("Sig10", 100.0), create("MISALLOC", 50.0));
    let sorted = super::tables::sort_signature_allocations(&[a2.clone(), a4.clone(), a3.clone(), a1.clone()]);
    assert_eq!(sorted, vec![a2, a1, a3, a4]);
}

#[test]
fn can_prefix_chromosomes() {
    // MiscUtilsTest
    assert_eq!(zero_prefixed("5"), "05");
    assert_eq!(zero_prefixed("5p13.2"), "05p13.2");
    assert_eq!(zero_prefixed("15"), "15");
    assert_eq!(zero_prefixed("15q21.1"), "15q21.1");
    assert_eq!(zero_prefixed("X"), "X");
}

#[test]
fn can_calc_fold_change() {
    // ExpressionsTest: formatFoldChange (cohort median), formatFoldChangeCancer, the field formats
    use super::tables::fold_change_display;
    assert_eq!(fold_change_display(1.0, Some(3.0)), "0.3");
    assert_eq!(fold_change_display(0.0, Some(0.0)), NOT_AVAILABLE);
    assert_eq!(fold_change_display(1.0, None), NOT_AVAILABLE);
    assert_eq!(tpm_field(1.234), "1.2");
    assert_eq!(two_digit(0.12345), "0.12");
    assert_eq!(two_digit(0.23445), "0.23");
}

#[test]
fn qc_status_interpretation() {
    // QcStatusInterpretationTest (the interpretations the report uses)
    use super::tables::{has_tumor_contaminated, is_fail_no_tumor};
    let qc = |s: &[&str]| crate::model::PurpleQC {
        status: s.iter().map(|x| x.to_string()).collect(),
        germline_aberrations: vec![],
        amber_mean_depth: 0,
        contamination: 0.0,
        total_copy_number_segments: 0,
        unsupported_copy_number_segments: 0,
        deleted_genes: 0,
    };
    assert!(has_tumor_contaminated(&qc(&["FAIL_CONTAMINATION"])));
    assert!(!has_tumor_contaminated(&qc(&["FAIL_NO_TUMOR"])));
    assert!(has_tumor_contaminated(&qc(&["PASS", "FAIL_CONTAMINATION"])));
    assert!(is_fail_no_tumor(&qc(&["FAIL_NO_TUMOR"])));
    assert!(!is_fail_no_tumor(&qc(&["FAIL_CONTAMINATION"])));
    assert!(is_fail_no_tumor(&qc(&["PASS", "FAIL_NO_TUMOR"])));
}

/// ReportWriterTest: reports from the Java tests' run (targeted, WGS tumor-only, WGS tumor/normal), in memory;
/// skipped without the Java test resources
#[test]
fn can_generate_test_reports_from_test_resources() {
    let Some(dir) = crate::tests::java_test_resources() else { return };
    let run = format!("{dir}/test_run");
    let virus = if std::path::Path::new(&format!("{run}/virus")).exists() { "virus" } else { "virusinterprtr" };
    let out = std::env::temp_dir().join(format!("orange_rs_report_test_{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    for (kind, reference) in [("PANEL", None), ("WGS", None), ("WGS", Some("ref_sample"))] {
        let mut args: Vec<String> = vec!["-experiment_type", kind, "-tumor", "tumor_sample", "-ref_genome_version", "37"].iter().map(|s| s.to_string()).collect();
        if let Some(r) = reference {
            args.extend(["-reference".to_string(), r.to_string(), "-linx_germline_dir".to_string(), format!("{run}/linx_germline")]);
        }
        for (opt, sub) in [
            ("-pipeline_version_file", "pipeline.version"),
            ("-purple_dir", "purple"),
            ("-purple_plot_dir", "purple/plot"),
            ("-linx_dir", "linx"),
            ("-linx_plot_dir", "linx/plot"),
            ("-lilac_dir", "lilac"),
            ("-chord_dir", "chord"),
            ("-cuppa_dir", "cuppa"),
            ("-peach_dir", "peach"),
            ("-sigs_dir", "sigs"),
            ("-virus_dir", virus),
            ("-isofox_dir", "isofox"),
        ] {
            args.extend([opt.to_string(), format!("{run}/{sub}")]);
        }
        args.extend(["-rna_sample".to_string(), "tumor_sample".to_string(), "-output_dir".to_string(), format!("{}/", out.display())]);
        let mut cb = hmf_common::config_builder::ConfigBuilder::new(crate::config::APP_NAME, crate::config::VERSION);
        crate::config::register(&mut cb);
        assert!(matches!(cb.parse(&args), Ok(hmf_common::config_builder::Parsed::Run)));
        let config = crate::config::OrangeConfig::from_builder(&cb).unwrap();
        crate::compat::with_java(false, || {
            let record = crate::algo::run(&config).unwrap();
            let doc = super::render_pdf(&record, &config).unwrap();
            assert!(doc.page_count() >= 5, "{kind} {reference:?}: {} pages", doc.page_count());
            assert!(doc.to_bytes().unwrap().starts_with(b"%PDF-"));
        });
    }
    std::fs::remove_dir_all(&out).ok();
}
