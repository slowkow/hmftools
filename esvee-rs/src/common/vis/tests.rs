//! Tests of the drawing and HTML helpers against Java. `testdata/expected.tsv` is the output of
//! `testdata/VisProbe.java`, which calls the hmf-common vis classes, JFreeSVG and j2html of the ESVEE jar
//! on the synthetic inputs repeated here (run on macOS, where "SansSerif" is Lucida Grande).

use std::collections::HashMap;

use super::html::*;
use super::svg::*;

fn expected() -> HashMap<&'static str, &'static str> {
    include_str!("testdata/expected.tsv").lines().filter_map(|line| line.split_once('\t')).collect()
}

fn check(name: &str, actual: &str) {
    let all = expected();
    let want = all.get(name).unwrap_or_else(|| panic!("no expected value {name}"));
    if *want != actual {
        let at = want.bytes().zip(actual.bytes()).position(|(a, b)| a != b).unwrap_or(want.len().min(actual.len()));
        let from = at.saturating_sub(60);
        panic!("{name} differs at byte {at}:\n java: {}\n rust: {}", &want[from..(at + 60).min(want.len())], &actual[from..(at + 60).min(actual.len())]);
    }
}

const RED: Color = Color(255, 0, 0);
const BLUE: Color = Color(0, 0, 255);

#[test]
fn css_builder_matches_java() {
    let border = CssBuilder::empty().border(CssSize::px(1.0), "solid", BLACK).border_collapse("collapse");
    let first_row = CssBuilder::empty()
        .border_left(CssSize::px(1.0), "solid", BLACK)
        .border_right(CssSize::px(1.0), "solid", BLACK)
        .border_top(CssSize::px(2.0), "solid", BLACK)
        .border_bottom(CssSize::px(1.0), "solid", BLACK)
        .border_collapse("collapse");
    check("css.empty", &format!("[{}]", CssBuilder::empty()));
    check("css.table", &border.margin_right(CssSize::px(10.0)).to_string());
    check("css.cell", &CssBuilder::empty().padding_left(CssSize::em(0.5)).padding_right(CssSize::em(0.5)).merge(&border).text_align("right").to_string());
    check(
        "css.left_header",
        &CssBuilder::empty()
            .background_color(LIGHT_GRAY)
            .writing_mode("vertical-rl")
            .transform("rotate(180deg)")
            .font_weight("bold")
            .text_align("center")
            .merge(&first_row)
            .to_string(),
    );
    check("css.merge_wins", &CssBuilder::empty().text_align("left").height(CssSize::em(1.0)).merge(&CssBuilder::empty().text_align("right")).to_string());
    check("css.div", &CssBuilder::empty().padding(CssSize::ZERO).margin(CssSize::ZERO).display("none").to_string());
    check("css.base_font", &base_font_style().to_string());
    check("css.sizes", &format!("{} {} {} {} {}", CssSize::px(0.00004), CssSize::px(0.00005), CssSize::em(1.23456), CssSize::px(-2.5), CssSize::ZERO));
    check("css.shadow", &CssBuilder::empty().left_box_shadow(CssSize::px(2.0), RED).right_box_shadow(CssSize::px(3.0), BLUE).font_size_pct(0.755).to_string());
}

#[test]
fn html_elements_match_j2html() {
    check("html.escape", &td("a<b>&\"c'd"));
    check("html.attrs", &tag("td", &[("rowspan", "3"), ("style", "color:#000000;")], &escape("x")));
    let rows = tag("tr", &[], &(td("") + &tag("td", &[], "<svg/>")));
    check("html.table", &styled_table(&rows, &CssBuilder::empty().border_spacing(CssSize::ZERO)));
    let js = javascript();
    // String.hashCode
    let hash = js.encode_utf16().fold(0i32, |h, c| h.wrapping_mul(31).wrapping_add(c as i32));
    let chars: Vec<char> = js.chars().collect();
    let head: String = chars[..40].iter().collect();
    let tail: String = chars[chars.len() - 20..].iter().collect();
    check("html.js", &format!("{} {} {} {}", js.encode_utf16().count(), hash, head.replace('\n', "\\n"), tail.replace('\n', "\\n")));
}

#[test]
fn read_info_table_matches_java() {
    let extra = [("Mismatch Info:", "1:2 <x>".to_string()), ("Trim count:", "0,3".to_string())];
    let full = ReadInfo {
        read_name: "read<1>",
        chromosome: "chrA",
        alignment_start: 100,
        alignment_end: 250,
        mate_alignment: Some(("chrB", 400, 0)),
        cigar: "10S140M",
        insert_size: Some(345),
        map_q: Some(60),
        extra_info: &extra,
        ..Default::default()
    };
    check("info.full", &render_read_info_table(&full));
    let unmapped = ReadInfo {
        read_name: "read2",
        chromosome: "chrA",
        alignment_start: 100,
        alignment_end: 250,
        cigar: "150M",
        mate_cigar: Some("100M"),
        orientation: Some("FR"),
        read_nm: Some(3),
        consensus_type: Some("DUAL"),
        consensus_read: Some("5;2;1"),
        second_map_q: Some(20),
        second_read_nm: Some(4),
        ..Default::default()
    };
    check("info.unmapped", &render_read_info_table(&unmapped));
}

#[test]
fn string_bounds_match_awt() {
    for text in ["9", "1,234,567", ">> chr7:1234567", "chrUn_1:99 <<", "INS", "A", "128"] {
        for size in [30, 26, 18] {
            let (width, height) = string_bounds(size, text);
            check(
                &format!("bounds.{size}.{text}"),
                &format!("0.0 {} {} {} {}", java_double(width), java_double(height), java_double(0.0 + width / 2.0), font_descent(size)),
            );
        }
    }
}

#[test]
fn java_double_to_string() {
    for (value, text) in [
        (0.0, "0.0"),
        (-0.0, "-0.0"),
        (12.0, "12.0"),
        (0.4, "0.4"),
        (0.001, "0.001"),
        (0.0001, "1.0E-4"),
        (1e7, "1.0E7"),
        (9999999.0, "9999999.0"),
        (12345678.9, "1.23456789E7"),
        (-0.39999999999999997, "-0.39999999999999997"),
        (1.0 / 3.0, "0.3333333333333333"),
        (f64::NAN, "NaN"),
        (f64::NEG_INFINITY, "-Infinity"),
    ] {
        assert_eq!(text, java_double(value));
    }
    assert_eq!("1,234,567", grouped(1234567));
    assert_eq!("999", grouped(999));
    assert_eq!("-1,000", grouped(-1000));
    assert_eq!("0", grouped(0));
}

fn cigar(text: &str) -> Vec<(i32, VisCigarOp)> {
    let mut out = Vec::new();
    let mut len = 0;
    for c in text.chars() {
        if let Some(d) = c.to_digit(10) {
            len = len * 10 + d as i32;
            continue;
        }
        let op = match c {
            'M' | '=' | 'X' => VisCigarOp::Match,
            'S' => VisCigarOp::S,
            'I' => VisCigarOp::I,
            'D' => VisCigarOp::D,
            'H' | 'N' => VisCigarOp::Skip,
            _ => VisCigarOp::Pad,
        };
        out.push((len, op));
        len = 0;
    }
    out
}

const READ_BASES: &[u8] = b"ACGTTGCANNACGTACGTAC";

fn read_quals() -> Vec<u8> {
    let mut quals: Vec<u8> = (0..READ_BASES.len()).map(|i| (i * 2 + 1) as u8).collect();
    quals[5] = 200; // negative as a Java byte
    quals
}

fn ref_model() -> BaseSeqViewModel {
    BaseSeqViewModel::from_bases(b"ACGTAGCATTACGGACGAACGTTT", 100)
}

fn read_model() -> BaseSeqViewModel {
    BaseSeqViewModel::create(101, &cigar("2S4M2I3M3D5M1H4S"), READ_BASES, &read_quals(), true).unwrap()
}

fn asm_model() -> BaseSeqViewModel {
    BaseSeqViewModel::from_string_with_cigar(b"ACGTAGCTTTACGG", &cigar("3S5M1D2I4M2N"), 98).unwrap()
}

fn opt(value: Option<bool>) -> String {
    value.map_or("null".to_string(), |v| v.to_string())
}

#[test]
fn base_seq_view_models_match_java() {
    let reference = ref_model();
    let mut read = read_model();
    check(
        "model.read",
        &format!(
            "{} {} {} {} {} {} {} {} {} {} {}",
            read.first_base_pos,
            read.last_base_pos,
            opt(read.left_is_forward_strand),
            opt(read.right_is_forward_strand),
            read.mismatch_count(&reference),
            read.get_base(101).soft_clip,
            read.get_base(106).right_insert_count,
            read.get_base(110).is_del(),
            read.get_base(118).is_missing(),
            read.get_base(99).is_missing(),
            read.get_base(104).base_q
        ),
    );
    // a quality above 127 is negative in Java
    assert_eq!(-56, read.get_base(106).base_q);
    read.clear_soft_clips();
    check("model.read_cleared", &format!("{} {}", read.mismatch_count(&reference), read.get_base(101).soft_clip));

    let asm = asm_model();
    check(
        "model.asm",
        &format!(
            "{} {} {} {} {} {} {}",
            asm.first_base_pos,
            asm.last_base_pos,
            asm.has_orientation(),
            asm.mismatch_count(&reference),
            asm.get_base(98).soft_clip,
            asm.get_base(106).is_del(),
            asm.get_base(105).right_insert_count
        ),
    );
    let empty = BaseSeqViewModel::from_bases(b"", 5);
    check("model.empty", &format!("{} {} {}", empty.first_base_pos, empty.last_base_pos, read.mismatch_count(&empty)));

    let exception = |r: super::VisResult<BaseSeqViewModel>| r.err().map_or("none".to_string(), |e| e.0.split(':').next().unwrap().to_string());
    check("model.overrun", &exception(BaseSeqViewModel::create(1, &cigar("30M"), READ_BASES, &read_quals(), false)));
    check("model.str_overrun", &exception(BaseSeqViewModel::from_string_with_cigar(b"ACGT", &cigar("5M"), 1)));
    let lead = BaseSeqViewModel::create(1, &cigar("3I2M"), b"ACGTA", &[1, 2, 3, 4, 5], false).unwrap();
    check("model.lead_insert", &format!("{} {} {}", lead.first_base_pos, lead.last_base_pos, lead.get_base(1).char_base().unwrap() as char));
    // qualities shorter than the bases fail as the array access does in Java
    assert!(BaseSeqViewModel::create(1, &cigar("5M"), b"ACGTA", &[1, 2, 3], false).is_err());
}

#[test]
fn base_sequences_render_as_jfreesvg() {
    let reference = ref_model();
    let read = read_model();
    let mut canvas = Canvas::new(12.0 * 40.0, 12.0);
    let init = canvas.transform();
    render_base_seq(&mut canvas, (0.0, 0.0), 12.0, BaseRegion::new(100, 112), &read, true, Some(&reference), true, false, false);
    canvas.set_transform(init);
    render_base_seq(&mut canvas, (15.0, 0.0), 12.0, BaseRegion::new(108, 123), &read, true, Some(&reference), false, true, true);
    check("svg.read", &canvas.svg_element());

    let asm = asm_model();
    let mut canvas = Canvas::new(12.0 * 30.0, 12.0);
    let init = canvas.transform();
    render_base_seq(&mut canvas, (0.0, 0.0), 12.0, BaseRegion::new(100, 110), &reference, false, None, true, true, false);
    canvas.set_transform(init);
    render_base_seq(&mut canvas, (13.0, 0.0), 12.0, BaseRegion::new(98, 112), &asm, false, Some(&reference), true, true, false);
    check("svg.ref_asm", &canvas.svg_element());

    let forward = BaseSeqViewModel::create(100, &cigar("5M"), b"ACGTT", &[37, 40, 0, 10, 20], false).unwrap();
    let mut canvas = Canvas::new(12.0 * 9.0, 12.0);
    render_base_seq(&mut canvas, (1.0, 0.0), 12.0, BaseRegion::new(100, 104), &forward, true, Some(&reference), true, true, false);
    check("svg.fwd", &canvas.svg_element());
    let mut canvas = Canvas::new(12.0 * 9.0, 12.0);
    render_base_seq(&mut canvas, (1.0, 0.0), 12.0, BaseRegion::new(100, 104), &forward, true, Some(&reference), true, true, true);
    check("svg.fwd_reversed_markers", &canvas.svg_element());

    let del = BaseSeqViewModel::create(100, &cigar("2M2I1M6D2M"), b"ACGGTAC", &[30; 7], false).unwrap();
    let mut canvas = Canvas::new(12.0 * 8.0, 12.0);
    render_base_seq(&mut canvas, (0.0, 0.0), 12.0, BaseRegion::new(102, 105), &del, false, Some(&reference), true, true, false);
    check("svg.del_end", &canvas.svg_element());
}

#[test]
fn coordinates_render_as_jfreesvg() {
    let scaling_factor = 12.0 / BASE_BOX_SIZE;
    let char_width = string_bounds(COORD_FONT_SIZE, "9").0;
    let max_string_width = string_bounds(COORD_FONT_SIZE, "1,000,012").0;
    let r1 = BaseRegion::new(999990, 1000012);
    let r2 = BaseRegion::new(2031, 2050);
    let boxes = r1.base_length() + 2 + r2.base_length() + 2;
    let mut canvas = Canvas::new(scaling_factor * BASE_BOX_SIZE * boxes as f64, scaling_factor * (max_string_width + 2.0 * char_width));
    let init = canvas.transform();
    let height = canvas.height;
    render_coords(&mut canvas, (0.0, 0.0), (scaling_factor * BASE_BOX_SIZE * (r1.base_length() + 2) as f64, height), 12.0, r1, 1000012, 10, max_string_width, false);
    canvas.set_transform(init);
    render_coords(
        &mut canvas,
        ((r1.base_length() + 2) as f64, 0.0),
        (scaling_factor * BASE_BOX_SIZE * (r2.base_length() + 2) as f64, height),
        12.0,
        r2,
        2031,
        10,
        max_string_width,
        true,
    );
    check("svg.coords", &canvas.svg_element());
}

#[test]
fn chromosome_labels_render_as_jfreesvg() {
    let label = |chromosome: Option<&str>, position, view: (i32, i32), reference: Option<(i32, i32)>, is_reversed, alignments: Vec<Alignment>| ChrLabel {
        chromosome: chromosome.map(|c| c.to_string()),
        position,
        view_region: BaseRegion::new(view.0, view.1),
        ref_region: reference.map(|r| BaseRegion::new(r.0, r.1)),
        is_reversed,
        alignments,
    };
    let labels = [
        label(Some("chr7"), 1500, (0, 24), Some((1001, 1500)), false, vec![Alignment::Right]),
        label(None, 0, (25, 31), None, false, vec![]),
        label(Some("12"), 88000, (32, 60), Some((88000, 88400)), true, vec![Alignment::Left]),
    ];
    check("svg.chr_labels", &render_chr_labels(12.0, &labels).unwrap().svg_element());
    let labels = [
        label(Some("chrX"), 20, (0, 30), Some((20, 300)), true, vec![Alignment::Right]),
        label(Some("chrX"), 900, (31, 70), Some((600, 900)), false, vec![Alignment::Left]),
    ];
    check("svg.chr_labels2", &render_chr_labels(12.0, &labels).unwrap().svg_element());
    assert!(render_chr_labels(12.0, &[]).is_err());
}

#[test]
fn colours_match_java() {
    let c = lighten(Color(151, 83, 9), 1.0 - 1.0 * 20.0 / 37.0);
    let o = overlapping_fragment_border_color();
    check("color.lighten", &format!("{},{},{} {},{},{}", c.0, c.1, c.2, o.0, o.1, o.2));
    assert_eq!(Color(151, 83, 9), lighten(Color(151, 83, 9), -0.5));
    assert_eq!(WHITE, lighten(Color(151, 83, 9), 1.5));
}

#[test]
fn base_region() {
    let region = BaseRegion::new(10, 19);
    assert_eq!(10, region.base_length());
    assert!(region.overlaps(&BaseRegion::new(19, 30)));
    assert!(region.overlaps(&BaseRegion::new(1, 10)));
    assert!(!region.overlaps(&BaseRegion::new(20, 30)));
    assert!(!region.overlaps(&BaseRegion::new(1, 9)));
}
