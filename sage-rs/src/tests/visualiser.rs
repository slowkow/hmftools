//! The visualiser. Java SAGE has no unit tests for its `vis` classes, so the expected values here were
//! produced by Java itself: `data/VisFixtures.java` calls the `hmf-common` and SAGE 5.1 classes on
//! synthetic input and prints one `key<TAB>value` line per case (`data/vis_java_fixtures.tsv.gz`).
//! The font-dependent values are those of Java on macOS.

use super::util::*;
use crate::types::*;
use crate::vis::*;
use crate::vis_aa::*;
use crate::vis_html::*;
use crate::vis_svg::*;
use std::collections::{BTreeSet, HashMap};
use std::io::Read as _;

fn fixtures() -> HashMap<String, String> {
    let mut text = String::new();
    flate2::read::GzDecoder::new(&include_bytes!("data/vis_java_fixtures.tsv.gz")[..]).read_to_string(&mut text).unwrap();
    text.lines().map(|l| l.split_once('\t').unwrap()).map(|(k, v)| (k.to_string(), v.to_string())).collect()
}

/// the synthetic chromosome of VisFixtures.java
fn genome(length: usize) -> Vec<u8> {
    let mut state: i64 = 12345;
    (0..length)
        .map(|_| {
            state = (state * 1103515245 + 12345) & 0x7fffffff;
            b"ACGT"[((state >> 16) & 3) as usize]
        })
        .collect()
}

#[test]
fn css_builder() {
    let f = fixtures();
    let e = CssBuilder::empty();
    assert_eq!(format!("[{e}]"), f["css.empty"]);
    let border = e.border(CssSize::px(1.0), "solid", BLACK).border_collapse("collapse");
    assert_eq!(border.margin_right(CssSize::px(10.0)).to_string(), f["css.table"]);
    let header = e.font_weight("bold").text_align("center").background_color(LIGHT_GRAY);
    let merged = e.padding_left(CssSize::em(0.5)).padding_right(CssSize::em(0.5)).merge(&border).merge(&header);
    assert_eq!(merged.to_string(), f["css.merge"]);
    assert_eq!(e.padding(CssSize::ZERO).margin(CssSize::ZERO).display("none").to_string(), f["css.hidden"]);
    let shadow = e.font_size_pt(8).height(CssSize::em(1.0)).left_box_shadow(CssSize::px(2.0), RED).right_box_shadow(CssSize::px(1.5), BLUE);
    assert_eq!(shadow.to_string(), f["css.shadow"]);
    assert_eq!(e.left_box_shadow(CssSize::px(2.0), RED).to_string(), f["css.shadow_only"]);
    assert_eq!(e.font_size_pct(0.123456).to_string(), f["css.pct"]);
    assert_eq!(e.text_align("left").color(PURPLE).merge(&e.text_align("right")).to_string(), f["css.override"]);
    let misc = e
        .border_bottom(CssSize::px(2.0), "solid", BLACK)
        .vertical_align("top")
        .writing_mode("vertical-rl")
        .border_spacing(CssSize::ZERO)
        .no_border()
        .font_family("sans-serif");
    assert_eq!(misc.to_string(), f["css.misc"]);
    // a builder is a value: deriving from it leaves it unchanged
    assert_eq!(border.to_string(), "border:1.0000px solid #000000;border-collapse:collapse;");
}

#[test]
fn css_size() {
    let f = fixtures();
    assert_eq!(CssSize::px(0.00004).to_string(), f["size.tiny"]);
    assert_eq!(CssSize::px(0.00005).to_string(), f["size.small"]);
    assert_eq!(CssSize::em(1.23456).to_string(), f["size.em"]);
    assert_eq!(CssSize::px(2.0).scale(-1.0).to_string(), f["size.neg"]);
    assert_eq!(CssSize::em(0.5f32 as f64).to_string(), f["size.half"]);
    assert!(CssSize::ZERO.is_zero());
}

#[test]
fn html_escaping_and_tags() {
    assert_eq!(escape("a<b>&\"c\"'d'"), "a&lt;b&gt;&amp;&quot;c&quot;&#x27;d&#x27;");
    assert_eq!(td("G>A"), "<td>G&gt;A</td>");
    assert_eq!(tag("td", &[("colspan", "9"), ("style", "a:b;")], "x"), "<td colspan=\"9\" style=\"a:b;\">x</td>");
    assert_eq!(styled_table("", &CssBuilder::empty()), "<table style=\"font-family:sans-serif;font-size:10pt;\"></table>");
}

#[test]
fn java_double_formatting() {
    for (value, expected) in [
        (0.4, "0.4"),
        (-0.0, "-0.0"),
        (0.0, "0.0"),
        (2136.0, "2136.0"),
        (1e-5, "1.0E-5"),
        (12345678.0, "1.2345678E7"),
        (1e7, "1.0E7"),
        (9999999.0, "9999999.0"),
        (0.001, "0.001"),
        (0.00099, "9.9E-4"),
        (253.49999999999997, "253.49999999999997"),
        (77.0 + 2.0 / 3.0, "77.66666666666667"),
        ((2.0f64 / 30.0) as f32 as f64, "0.06666667014360428"),
        (-75.884765625, "-75.884765625"),
    ] {
        assert_eq!(java_double(value), expected);
    }
    assert_eq!(grouped(7920552), "7,920,552");
    assert_eq!(grouped(999), "999");
    assert_eq!(grouped(1000), "1,000");
}

#[test]
fn font_metrics() {
    let f = fixtures();
    let mut checked = 0;
    for size in [18, 26, 30] {
        for text in ["7,920,552", "G", "12", "*", "W", ""] {
            let (width, height) = string_bounds(size, text);
            let actual = format!("{} {} {}", java_double(width), java_double(height), font_descent(size));
            assert_eq!(actual, f[&format!("font.{size}.{text}")], "size {size} text {text}");
            checked += 1;
        }
    }
    assert_eq!(checked, 18);
}

#[test]
fn colours() {
    let f = fixtures();
    let show = |c: Color| format!("{},{},{}", c.0, c.1, c.2);
    assert_eq!(show(interpolate_colors(RED, GREEN, 37.0 / 60.0)), f["color.interpolate"]);
    assert_eq!(show(lighten(Color(151, 83, 9), 1.0 - 1.0 * 11.0 / 37.0)), f["color.lighten"]);
    assert_eq!(show(overlapping_fragment_border_color()), f["color.overlap"]);
}

#[test]
fn svg_box_and_coordinates() {
    let f = fixtures();
    assert_eq!(render_colored_box(12.0, Color(1, 2, 3)).svg_element(), f["svg.box"]);
    assert_eq!(render_coords(12.0, BaseRegion::new(995, 1012), 1003, 10).svg_element(), f["svg.coords"]);
}

fn test_borders() -> Vec<(i32, Vec<(BorderLocation, Color)>)> {
    vec![
        (115, vec![(BorderLocation::Left, BLACK), (BorderLocation::Right, BLACK)]),
        (112, vec![(BorderLocation::Left, BLUE)]),
        (118, vec![(BorderLocation::Right, BLUE)]),
        (104, vec![(BorderLocation::Left, DARK_BLUE)]),
        (126, vec![(BorderLocation::Right, DARK_BLUE)]),
    ]
}

#[test]
fn svg_base_sequences() {
    let f = fixtures();
    let reference = genome(2000);
    let view = BaseRegion::new(100, 130);
    let ref_model = BaseSeqViewModel::from_bases(&reference[89..140], 90);
    let borders = test_borders();
    assert_eq!(render_base_seq(12.0, view, &ref_model, false, &borders, None).svg_element(), f["svg.ref"]);

    // soft clips, an insert, a delete, mismatches and shaded qualities on the reverse strand
    let bases = &f["read1.bases"];
    let quals: Vec<u8> = (0..bases.len()).map(|i| (5 + (i * 7) % 36) as u8).collect();
    let mut read = build_sam_record_quals(102, "3S8M2I6M3D9M2S", bases, quals);
    read.flags |= FLAG_REVERSE;
    assert_eq!(unclipped_start(&read), 99);
    assert_eq!(unclipped_end(&read), 129);
    let read_model = BaseSeqViewModel::from_read(&read);
    assert_eq!(render_base_seq(12.0, view, &read_model, true, &borders, Some(&ref_model)).svg_element(), f["svg.read1"]);

    // two overlapping mates and the read formed from them
    let str_of = |b: &[u8]| String::from_utf8(b.to_vec()).unwrap();
    let first = build_sam_record_quals(96, "20M", &str_of(&reference[95..115]), vec![30; 20]);
    let mut second = build_sam_record_quals(113, "5S18M", &format!("CCCCC{}", str_of(&reference[112..130])), vec![30; 23]);
    second.flags |= FLAG_REVERSE;
    let consensus = build_sam_record_quals(96, "35M", &str_of(&reference[95..130]), vec![37; 35]);
    let first_model = BaseSeqViewModel::from_read(&first);
    let second_model = BaseSeqViewModel::from_read(&second);
    let fragment_model = BaseSeqViewModel::from_consensus_fragment(&consensus, Some(&first_model), Some(&second_model));
    assert_eq!(fragment_model.left_is_forward_strand, Some(true));
    assert_eq!(fragment_model.right_is_forward_strand, Some(false));
    assert_eq!(render_base_seq(12.0, view, &fragment_model, true, &borders, Some(&ref_model)).svg_element(), f["svg.fragment"]);
    assert_eq!(render_base_seq(12.0, view, &first_model, true, &borders, Some(&ref_model)).svg_element(), f["svg.first"]);
    assert_eq!(render_base_seq(12.0, view, &second_model, true, &borders, Some(&ref_model)).svg_element(), f["svg.second"]);
}

#[test]
fn svg_gene_tracks() {
    let f = fixtures();
    let mut insert_before = GeneRegionViewModel::amino_acid(BaseRegion::new(114, 116), 2, b'W', b'K', false);
    insert_before.bases_before_right_insert = 2;
    let models = vec![
        GeneRegionViewModel::intronic(BaseRegion::new(90, 104)),
        GeneRegionViewModel::non_coding(BaseRegion::new(105, 110)),
        GeneRegionViewModel::amino_acid(BaseRegion::new(111, 113), 1, b'M', b'M', true),
        insert_before,
        GeneRegionViewModel::del(BaseRegion::new(117, 119)),
        GeneRegionViewModel::amino_acid(BaseRegion::new(120, 122), 3, b'Q', b'*', false),
        GeneRegionViewModel::amino_acid(BaseRegion::new(123, 125), 3, b'I', b'I', false),
    ];
    for pos_strand in [true, false] {
        let (ref_canvas, alt_canvas) = render_gene_data(12.0, BaseRegion::new(100, 124), pos_strand, &models, &models);
        assert_eq!(ref_canvas.svg_element(), f[&format!("svg.gene.{pos_strand}.ref")]);
        assert_eq!(alt_canvas.svg_element(), f[&format!("svg.gene.{pos_strand}.alt")]);
    }
}

fn show_events(events: &[AminoAcidEvent]) -> String {
    events
        .iter()
        .map(|e| match *e {
            AminoAcidEvent::Substitution { pos, alt } => format!("S:{}:{}", pos, alt as char),
            AminoAcidEvent::Ins { pos, length } => format!("I:{pos}:{length}"),
            AminoAcidEvent::Del { start, end } => format!("D:{start}:{end}"),
            AminoAcidEvent::Ext { pos } => format!("E:{pos}"),
            AminoAcidEvent::StartLost { pos } => format!("L:{pos}"),
        })
        .collect::<Vec<_>>()
        .join(" ")
}

#[test]
fn amino_acid_event_parsing() {
    let f = fixtures();
    let mut checked = 0;
    for (key, expected) in f.iter().filter(|(k, _)| k.starts_with("hgvs.")) {
        let hgvs = &key["hgvs.".len()..];
        let actual = match parse_amino_acid_events(hgvs) {
            Ok(events) => format!("ok {}", show_events(&events)),
            Err(_) => "err".to_string(),
        };
        assert_eq!(&actual, expected, "{hgvs}");
        checked += 1;
    }
    assert_eq!(checked, 45);
    assert_eq!(parse_amino_acid_events("p.Gly12_Ala14delinsVal").unwrap(), vec![
        AminoAcidEvent::Del { start: 12, end: 13 },
        AminoAcidEvent::Substitution { pos: 14, alt: b'V' }
    ]);
}

fn show_models(models: &[GeneRegionViewModel]) -> String {
    models
        .iter()
        .map(|m| {
            let (s, e) = (m.region.start, m.region.end);
            match m.kind {
                GeneRegionKind::AminoAcid { amino_acid_pos, ref_, alt, is_start } => {
                    format!("A:{s}-{e}:{}:{}:{}:{}:{}", amino_acid_pos, ref_ as char, alt as char, m.bases_before_right_insert, is_start)
                }
                GeneRegionKind::NonCodingExonic => format!("N:{s}-{e}:{}", m.bases_before_right_insert),
                GeneRegionKind::Intronic => format!("I:{s}-{e}:{}", m.bases_before_right_insert),
                GeneRegionKind::Del => format!("D:{s}-{e}"),
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

const AMINO_ACIDS: &str = "MACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVW";

fn test_transcript(name: &str) -> Transcript {
    let (pos_strand, coding_start, coding_end) = match name {
        "pos" => (true, 1110, 1280),
        "neg" => (false, 1110, 1280),
        "posfull" => (true, 1100, 1300),
        "neg2" => (false, 1161, 1245),
        _ => panic!("{name}"),
    };
    Transcript {
        gene_id: "ENSG0001".to_string(),
        gene_name: "GENE1".to_string(),
        trans_name: "ENST0001".to_string(),
        pos_strand,
        coding_start: Some(coding_start),
        coding_end: Some(coding_end),
        exons: vec![(1100, 1130, -1), (1160, 1210, -1), (1240, 1300, -1)],
        amino_acids: AMINO_ACIDS.as_bytes().to_vec(),
    }
}

#[test]
fn gene_regions_and_labels() {
    let f = fixtures();
    for name in ["pos", "neg", "posfull", "neg2"] {
        let transcript = test_transcript(name);
        let regions = gene_regions(&transcript).unwrap();
        assert_eq!(show_models(&regions), f[&format!("regions.{name}")], "{name}");
        assert_eq!(show_models(&ref_gene_region_view_models(&transcript, &regions)), f[&format!("refmodels.{name}")], "{name}");
        // amino acid positions beyond the protein are left out
        let mut short = transcript.clone();
        short.amino_acids.truncate(30);
        assert_eq!(show_models(&ref_gene_region_view_models(&short, &regions)), f[&format!("refmodels.short.{name}")], "{name}");
        let labels: Vec<String> = [1000, 1100, 1115, 1130, 1131, 1159, 1160, 1210, 1220, 1240, 1300, 1301, 1400]
            .iter()
            .map(|p| gene_region_label(&transcript, *p))
            .collect();
        assert_eq!(labels.join("|"), f[&format!("labels.{name}")], "{name}");
    }
    let mut non_coding = test_transcript("pos");
    non_coding.coding_start = None;
    assert!(gene_regions(&non_coding).is_err());
}

/// Every combination of four transcripts, six variants, eighteen protein changes and two view regions
#[test]
fn alt_gene_region_view_models_match_java() {
    let f = fixtures();
    let reference = genome(2000);
    let ref_bases = |start: i32, end: i32| -> Result<Vec<u8>, String> {
        if start < 1 || end < start - 1 || end as usize > reference.len() {
            return Err(format!("invalid range {start}-{end}"));
        }
        Ok(reference[(start - 1) as usize..end as usize].to_vec())
    };
    let mut checked = 0;
    let mut failures = 0;
    let mut svgs = 0;
    for (key, expected) in f.iter().filter(|(k, _)| k.starts_with("altmodels.")) {
        // altmodels.<transcript>.<pos><ref>><alt>.<hgvs>.<view start>
        let rest = &key["altmodels.".len()..];
        let (name, rest) = rest.split_once('.').unwrap();
        let (variant_text, rest) = rest.split_once('.').unwrap();
        let (hgvs, view_start) = rest.rsplit_once('.').unwrap();
        let (position_and_ref, alt) = variant_text.split_once('>').unwrap();
        let digits = position_and_ref.chars().take_while(|c| c.is_ascii_digit()).count();
        let variant = SimpleVariant::new(position_and_ref[..digits].parse().unwrap(), position_and_ref[digits..].as_bytes(), alt.as_bytes());
        let view = if view_start == "1105" { BaseRegion::new(1105, 1281) } else { BaseRegion::new(1170, 1200) };
        let transcript = test_transcript(name);
        let regions = gene_regions(&transcript).unwrap();
        let alt_context = AltContext { variant: &variant, chromosome_length: 2000, ref_bases: &ref_bases };
        let models = parse_amino_acid_events(hgvs)
            .and_then(|events| alt_gene_region_view_models(&transcript, &regions, &events, view, &alt_context));
        let actual = match &models {
            Ok(models) => format!("ok {}", show_models(models)),
            Err(_) => "err".to_string(),
        };
        assert_eq!(&actual, expected, "{key}");
        checked += 1;
        failures += models.is_err() as usize;
        if let (Some(svg), Ok(models)) = (f.get(&format!("svg.{key}")), &models) {
            let ref_models = ref_gene_region_view_models(&transcript, &regions);
            let (ref_canvas, alt_canvas) = render_gene_data(12.0, view, transcript.pos_strand, &ref_models, models);
            assert_eq!(&(ref_canvas.svg_element() + &alt_canvas.svg_element()), svg, "{key}");
            svgs += 1;
        }
    }
    assert_eq!(checked, 4 * 6 * 18 * 2);
    assert_eq!(failures, 132);
    assert_eq!(svgs, 60);
}

#[test]
fn read_info_table() {
    let f = fixtures();
    let names = vec!["chr1".to_string(), "chr2".to_string()];
    let mut first = build_sam_record(1000, "100M", &"A".repeat(100));
    first.name = "read<1>&\"x\"'".to_string();
    first.flags = FLAG_PAIRED | FLAG_PROPER | FLAG_FIRST | FLAG_MATE_REVERSE;
    first.mate_tid = 1;
    first.mate_start = 1200;
    first.mate_cigar = Some("5S95M".to_string());
    first.insert_size = -295;
    first.nm = Some(3);
    first.has_consensus_tag = true;
    first.consensus = ConsensusType::Dual;
    first.consensus_read_tag = Some("7;3;4".to_string());
    let mut second = build_sam_record(1200, "5S95M", &"A".repeat(100));
    second.mapq = 59;
    second.nm = None;
    assert_eq!(render_read_info_table(&first, Some(&second), "chr1", &names), f["readinfo.paired"]);

    let mut single = create_sam_record_unpaired("r2", 1000, &"A".repeat(60), "10S50M");
    single.mapq = 13;
    single.nm = None;
    assert_eq!(render_read_info_table(&single, None, "chr1", &names), f["readinfo.single"]);

    // no mate cigar: the mate's end is unknown
    let mut no_mate_cigar = build_sam_record(1000, "50M", &"A".repeat(50));
    no_mate_cigar.name = "r3".to_string();
    no_mate_cigar.flags = FLAG_PAIRED | FLAG_FIRST | FLAG_REVERSE;
    no_mate_cigar.mate_tid = 0;
    no_mate_cigar.mate_start = 900;
    no_mate_cigar.insert_size = -150;
    no_mate_cigar.mapq = 13;
    no_mate_cigar.nm = Some(2);
    no_mate_cigar.has_consensus_tag = true;
    assert_eq!(render_read_info_table(&no_mate_cigar, None, "chr1", &names), f["readinfo.nomc"]);
}

#[test]
fn filter_names_follow_java_hash_set_order() {
    let f = fixtures();
    let mut checked = 0;
    for (key, expected) in f.iter().filter(|(k, _)| k.starts_with("filters.")) {
        let names: Vec<&str> = key["filters.".len()..].split(',').collect();
        assert_eq!(&filter_names_in_java_set_order(&names).join(","), expected);
        // the order does not depend on the order the names are given in, unless two share a bucket
        let mut reversed = names.clone();
        reversed.reverse();
        let reversed_order = filter_names_in_java_set_order(&reversed).join(",");
        if key.contains("dedupSnvMnv") {
            assert_eq!(reversed_order, "dedupIndel,dedupSnvMnv");
        } else {
            assert_eq!(&reversed_order, expected);
        }
        checked += 1;
    }
    assert_eq!(checked, 4);
}

#[test]
fn file_names() {
    let variant = SimpleVariant::new(100, b"A", b"T");
    let one: BTreeSet<String> = ["5_3_7_3_ACGTACGTACG".to_string()].into();
    assert_eq!(form_filename("chr1", &variant, "TUMOR", None, None, &one, "5_3_7_3_ACGTACGTACG"), "TUMORchr1_100_A_T.html.gz");
    assert_eq!(
        form_filename("chr1", &variant, "TUMOR", Some("KRAS"), Some("missense_variant"), &one, "5_3_7_3_ACGTACGTACG"),
        "TUMOR.sage.KRAS_chr1_100_A_T_missense.html.gz"
    );
    // the gene form always has the chr prefix; the plain form keeps the chromosome as given
    assert_eq!(
        form_filename("1", &variant, "TUMOR", Some("KRAS"), Some("splice_donor_variant&intron_variant"), &one, "x"),
        "TUMOR.sage.KRAS_chr1_100_A_T_splice_donor_intron.html.gz"
    );
    assert_eq!(form_filename("1", &variant, "TUMOR", None, None, &one, "x"), "TUMOR1_100_A_T.html.gz");
    assert_eq!(
        form_filename("chr1", &variant, "T", Some("G"), Some("stop gained & other_variant"), &one, "x"),
        "T.sage.G_chr1_100_A_T_stopgained_other.html.gz"
    );
    // several read contexts for one variant: numbered in the order of their keys
    let three: BTreeSet<String> = ["a", "b", "c"].iter().map(|s| s.to_string()).collect();
    assert_eq!(form_filename("chr1", &variant, "TUMOR", None, None, &three, "a"), "TUMORchr1_100_A_T_1.html.gz");
    assert_eq!(form_filename("chr1", &variant, "TUMOR", None, Some("x"), &three, "c"), "TUMORchr1_100_A_T_3_x.html.gz");
    let ten: BTreeSet<String> = (0..10).map(|i| format!("k{i}")).collect();
    assert_eq!(form_filename("chr1", &variant, "TUMOR", None, None, &ten, "k1"), "TUMORchr1_100_A_T_02.html.gz");
    assert_eq!(variant_key("chr1", &variant), "chr1_100_A_T");
}

#[test]
fn downsample_counts() {
    // ReadContextMatch ordinals: NONE, CORE, PARTIAL_CORE, FULL, REALIGNED, REF, SIMPLE_ALT
    let defaults: Vec<usize> = (0..7).map(|t| max_reads_shown(t, 0)).collect();
    assert_eq!(defaults, vec![5, 10, 10, 40, 20, 20, 5]);
    assert!((0..7).all(|t| max_reads_shown(t, -1) == 1000));
    assert!((0..7).all(|t| max_reads_shown(t, 7) == 7));
    assert!((0..7).all(|t| max_reads_shown(t, 5000) == 1000));
    for (count, max_reads) in [(0, 5), (3, 5), (5, 5), (6, 5), (76, 40), (397, 5), (1001, 1000), (5000, 7)] {
        let indices = downsample_indices(count, max_reads);
        assert_eq!(indices.len(), count.min(max_reads));
        assert!(indices.windows(2).all(|w| w[0] < w[1]), "distinct and in arrival order");
        assert!(indices.iter().all(|i| *i < count));
        // the same every time
        assert_eq!(indices, downsample_indices(count, max_reads));
    }
    assert_eq!(downsample_indices(10, 5), vec![0, 2, 4, 6, 8]);
}

#[test]
fn vis_variants_option() {
    let variants = parse_vis_variants("chr21:5234885:G:A;chr21:7920622:G:GTCCAT").unwrap();
    assert_eq!(variants.len(), 2);
    assert_eq!(variants[1], ("chr21".to_string(), SimpleVariant::new(7920622, b"G", b"GTCCAT")));
    assert!(parse_vis_variants("chr21:5234885:G").is_err());
    assert!(parse_vis_variants("chr21:x:G:A").is_err());
    assert!(parse_vis_variants("chr21:5:G:A;").is_err());

    let mut cfg = VisConfig::default();
    assert!(!cfg.enabled());
    assert!(cfg.process_variant("chr1", &variants[0].1));
    cfg.specific_variants = variants.clone();
    assert!(cfg.enabled() && cfg.has_variants() && !cfg.visualiser_only_mode());
    assert!(cfg.process_variant("chr21", &variants[0].1));
    assert!(!cfg.process_variant("chr20", &variants[0].1));
    assert!(!cfg.process_variant("chr21", &SimpleVariant::new(5234885, b"G", b"T")));
    let pass_only = VisConfig { pass_only: true, ..Default::default() };
    assert!(pass_only.enabled() && !pass_only.has_variants());
}

#[test]
fn each_vis_variant_adds_a_region() {
    let dir = std::env::temp_dir().join(format!("sage_rs_vis_cfg_{}", std::process::id()));
    let vcf = dir.join("out.sage.vcf.gz");
    let mut args: Vec<String> = ["-tumor", "T", "-ref_genome", "none.fasta", "-ref_genome_version", "38", "-output_vcf"].iter().map(|s| s.to_string()).collect();
    args.push(vcf.to_string_lossy().to_string());
    args.extend(["-vis_variants", "chr21:1000:G:A;chr21:1300:C:T;chr5:500:A:AT", "-vis_max_support_reads", "-1"].iter().map(|s| s.to_string()));
    let cfg = crate::config::Config::from_args(&args).unwrap();
    // the second variant's region overlaps the first, which is widened
    assert_eq!(cfg.specific_regions, vec![("chr21".to_string(), 800, 1500), ("chr5".to_string(), 300, 700)]);
    assert_eq!(cfg.specific_chr, vec!["chr21".to_string(), "chr5".to_string()]);
    assert_eq!(cfg.vis.max_support_reads, -1);
    assert_eq!(cfg.vis.output_dir, dir.join("vis").to_string_lossy());
    assert!(dir.join("vis").is_dir());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn reportable_variants_of_a_purple_vcf() {
    let path = std::env::temp_dir().join(format!("sage_rs_vis_purple_{}.vcf", std::process::id()));
    let impact = "GENE1,ENST0001,missense_variant,MISSENSE,false,c.35G>T,p.Gly12Val,,MISSENSE,1";
    std::fs::write(
        &path,
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT\n\
             chr1\t100\t.\tG\tT\t100\tPASS\tTIER=PANEL;IMPACT={impact};REPORTED\tGT\t0/1\n\
             chr1\t200\t.\tG\tT\t100\tPASS\tIMPACT={impact}\tGT\t0/1\n\
             chr2\t300\t.\tGA\tG,GAA\t100\tPASS\tREPORTED;IMPACT=a,b\tGT\t0/1\n"
        ),
    )
    .unwrap();
    let reportable = load_reportable_variants(&path.to_string_lossy()).unwrap();
    std::fs::remove_file(&path).unwrap();
    assert_eq!(reportable.len(), 2);
    assert_eq!((reportable[0].chromosome.as_str(), &reportable[0].variant), ("chr1", &SimpleVariant::new(100, b"G", b"T")));
    assert_eq!(reportable[1].variant, SimpleVariant::new(300, b"GA", b"G"));
    let parsed = impact_from_values(&reportable[0].impact_values).unwrap();
    assert_eq!(
        parsed,
        VariantImpact {
            gene_name: "GENE1".to_string(),
            canonical_transcript: "ENST0001".to_string(),
            canonical_effect: "missense_variant".to_string(),
            canonical_hgvs_protein: "p.Gly12Val".to_string(),
        }
    );
    // anything but ten values gives an empty impact, as VariantImpactSerialiser does
    assert_eq!(impact_from_values(&reportable[1].impact_values).unwrap(), VariantImpact::default());
    let mut bad = reportable[0].impact_values.clone();
    bad[3] = "MISSENSE_OR_SO".to_string();
    assert!(impact_from_values(&bad).is_err());
}

#[test]
fn context_view_models() {
    let ctx = |variant: SimpleVariant, read_bases: &str, var_index: i32| crate::read_context::VariantReadContext {
        variant,
        alignment_start: 0,
        alignment_end: 0,
        ref_bases: Vec::new(),
        read_bases: read_bases.as_bytes().to_vec(),
        core_index_start: var_index - 2,
        var_index,
        core_index_end: var_index + 2,
        homology: None,
        max_repeat: None,
        all_repeats: Vec::new(),
        core_position_start: 0,
        core_position_end: 0,
        read_cigar: String::new(),
        has_indel_in_core: false,
        artefact: None,
        max_ref_repeat: None,
        extended_ref_bases: None,
        marked_invalid: false,
    };
    let kinds = |m: &BaseSeqViewModel| -> String {
        m.bases
            .iter()
            .map(|b| match b.kind {
                BaseKind::Char(c) => format!("{}{}", c as char, if b.right_insert_count > 0 { b.right_insert_count.to_string() } else { String::new() }),
                BaseKind::Del => "-".to_string(),
                BaseKind::Missing => "?".to_string(),
            })
            .collect()
    };
    let snv = context_view_model_from_variant(&ctx(SimpleVariant::new(100, b"A", b"T"), "ACGTTGCA", 4));
    assert_eq!((snv.pos_start, kinds(&snv)), (96, "ACGTTGCA".to_string()));
    let del = context_view_model_from_variant(&ctx(SimpleVariant::new(100, b"TGG", b"T"), "ACGTTGCA", 4));
    assert_eq!((del.pos_start, kinds(&del)), (96, "ACGTT--GCA".to_string()));
    let ins = context_view_model_from_variant(&ctx(SimpleVariant::new(100, b"T", b"TGC"), "ACGTTGCA", 4));
    assert_eq!((ins.pos_start, kinds(&ins)), (96, "ACGTT2A".to_string()));
    assert_eq!(indexed_bases_key(&ctx(SimpleVariant::new(100, b"A", b"T"), "ACGTTGCA", 4)), "4_2_6_2_ACGTTGCA");
}
