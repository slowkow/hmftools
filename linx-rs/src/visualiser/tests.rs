//! Tests of the visualiser's data steps (synthetic data).

use super::circos_data::*;
use super::data::*;
use super::genome::*;

#[test]
fn shorthand_as_java() {
    assert_eq!(shorthand(99), "99");
    assert_eq!(shorthand(1500), "1.5k");
    assert_eq!(shorthand(99_949), "99.9k");
    assert_eq!(shorthand(99_950), "0.1m");
    assert_eq!(shorthand(2_345_678), "2.3m");
}

#[test]
fn log_distance_and_scale_contig() {
    assert_eq!(log_distance(1), 10);
    assert_eq!(log_distance(1000), 37);
    assert_eq!(log_distance(1_000_000), 226);
    let c = ScaleContig::new("1", &[100, 1100, 1_001_100]);
    assert_eq!(c.scale(100).unwrap(), 1);
    assert_eq!(c.scale(1100).unwrap(), 38);
    assert_eq!(c.scale(1_001_100).unwrap(), 38 + 226);
    assert_eq!(c.length(), 264);
    assert!(c.scale(5).is_err());
    // interpolation between keys, floored; clamped outside them
    assert_eq!(c.interpolate(600), 1 + (0.5f64 * 37.0).floor() as i32);
    assert_eq!(c.interpolate(50), 1);
    assert_eq!(c.interpolate(2_000_000), 264);
}

#[test]
fn cna_axis_positions_as_java() {
    assert_eq!(cna_axis_positions(2), "0.5r,1.0r");
    assert_eq!(cna_axis_positions(3), "0.3333r,0.6667r,1.0r");
    let p = cna_axis_positions(20);
    assert_eq!(p.split(',').count(), 7 + 2);
    assert!(p.ends_with("0.4r,0.9r"));
}

#[test]
fn protein_domain_hues_as_java() {
    // Color.getHSBColor(hue(i, 3), 0.5, 0.8392) for three domains
    let c: Vec<(u8, u8, u8)> = (0..3).map(|i| hsb_to_rgb(protein_hue(i, 3), 0.5, 0.8392)).collect();
    assert_eq!(protein_hue(0, 3), 21f32 / 360f32 + 1f32 / 3f32 / 2f32);
    assert_eq!(c.len(), 3);
    assert_eq!(hsb_to_rgb(0.0, 0.0, 0.8392), (214, 214, 214));
    assert_eq!(hsb_to_rgb(0.0, 1.0, 1.0), (255, 0, 0));
}

#[test]
fn downsample_targets() {
    assert_eq!(group_downsample_targets(&[0, 100, 300], 50, 5), vec![0, 5 + 10, 5 + 30]);
    let l: Vec<i32> = (0..10).collect();
    assert_eq!(downsample_list(&l, 4), vec![0, 2, 5, 7]);
    let r = regions_from_positions(&[GenomePosition::new("chr1", 100), GenomePosition::new("chr1", 50)]);
    assert_eq!(r, vec![Region::new("chr1", 0, 50), Region::new("chr1", 50, 100), Region::new("chr1", 100, i32::MAX)]);
}

#[test]
fn cobalt_centering() {
    let pos = vec![GenomePosition::new("1", 1200), GenomePosition::new("1", 1400), GenomePosition::new("1", 1900)];
    let r = center_cobalt_positions(&[("1".into(), 1001, 1.0), ("1".into(), 5001, 1.0), ("1".into(), 1, 1.0)], &pos);
    assert_eq!(r[0].1, 1400); // three breakpoints: the middle one
    assert_eq!(r[1].1, 5501); // none: the window's middle
    let r = center_cobalt_positions(&[("1".into(), 1001, 1.0)], &pos[..2]);
    assert_eq!(r[0].1, 1300); // two: the middle of the segment between them
}

#[test]
fn chromosome_order_and_names() {
    let mut p = vec![GenomePosition::new("chrX", 5), GenomePosition::new("chr10", 3), GenomePosition::new("chr2", 9), GenomePosition::new("chr2", 1)];
    sort_positions(&mut p);
    let names: Vec<String> = p.iter().map(|x| format!("{}:{}", x.chromosome, x.position)).collect();
    assert_eq!(names, ["chr2:1", "chr2:9", "chr10:3", "chrX:5"]);
    assert_eq!(circos_contig("chr9"), "hs9");
    assert!(human_chromosome_contains("chrY"));
    assert!(!human_chromosome_contains("chrM"));
    assert!(!human_chromosome_contains("23"));
}

#[test]
fn legend_order_as_r_sorts() {
    let mut v = vec!["ANK_REPEAT", "ABD", "3-Prime UTR", "ANK_REP_REGION", "abc", "ZINC_FINGER_C2H2_2"];
    v.sort_by(|a, b| super::fusion::r_collate(a, b));
    // sort(c(...)) in the pipeline's R (ICU collation)
    assert_eq!(v, ["3-Prime UTR", "abc", "ABD", "ANK_REP_REGION", "ANK_REPEAT", "ZINC_FINGER_C2H2_2"]);
}

#[test]
fn merge_count_order_as_java() {
    // 13 keys: still 16 buckets after the 13th insertion (merge resizes at the start of the next call)
    let keys: Vec<i32> = (1..=13).map(|i| i * 16 + 1).collect(); // all in bucket 1 of 16; buckets 1 and 17 of 32
    let o: Vec<i32> = merge_count_order(&keys).iter().map(|x| x.0).collect();
    let expect: Vec<i32> = keys.iter().rev().copied().collect();
    assert_eq!(o, expect);
    // one more call (an existing key) resizes to 32 buckets
    let mut more = keys.clone();
    more.push(17);
    let o: Vec<(i32, i64)> = merge_count_order(&more);
    let even: Vec<i32> = keys.iter().filter(|k| (**k & 31) == 1).rev().copied().collect();
    assert_eq!(o.iter().map(|x| x.0).take(even.len()).collect::<Vec<_>>(), even);
    assert_eq!(o.iter().find(|x| x.0 == 17).unwrap().1, 2);
}

// ---------------------------------------------------------------------------------------------------------
// LINX's Java visualiser tests, ported

use super::files::*;

#[test]
fn circos_config_writer_cna_axis_positions() {
    // CircosConfigWriterTest.testCnaAxisPositions
    let max_tracks = 42;
    let s = cna_axis_positions(max_tracks);
    let cn: Vec<String> = s.split(',').map(|p| (java_round(max_tracks as f64 * p.replace('r', "").parse::<f64>().unwrap()) + 2).to_string()).collect();
    assert_eq!(cn.join(","), "3,4,5,6,7,8,9,10,20,30,40");
}

#[test]
fn circos_data_writer_shorthand() {
    // CircosDataWriterTest.testShorthand
    assert_eq!(shorthand(1), "1");
    assert_eq!(shorthand(99), "99");
    assert_eq!(shorthand(100), "0.1k");
    assert_eq!(shorthand(1_000), "1.0k");
    assert_eq!(shorthand(99_949), "99.9k");
    assert_eq!(shorthand(99_950), "0.1m");
    assert_eq!(shorthand(100_000), "0.1m");
}

#[test]
fn colour_picker_contig_colours() {
    // ColorPickerTest.testContigColours
    assert_ne!(hex_contig_colour("1"), "#000000");
    assert_eq!(hex_contig_colour("1"), hex_contig_colour("chr1"));
    assert_eq!(hex_contig_colour("MT"), "#000000");
}

#[test]
fn protein_domain_colours_hue() {
    // ProteinDomainColorsTest.testHue
    let close = |a: f32, b: f32| assert!((a - b).abs() < 0.1, "{a} {b}");
    close(21.0 + 90.0, 360.0 * protein_hue(0, 1));
    close(21.0 + 90.0, 360.0 * protein_hue(0, 2));
    close(21.0 + 270.0, 360.0 * protein_hue(1, 2));
    close(21.0 + 60.0, 360.0 * protein_hue(0, 3));
    close(21.0 + 180.0 + 60.0, 360.0 * protein_hue(1, 3));
    close(21.0 + 120.0, 360.0 * protein_hue(2, 3));
    close(21.0 + 60.0, 360.0 * protein_hue(0, 4));
    close(21.0 + 180.0 + 60.0, 360.0 * protein_hue(1, 4));
    close(21.0 + 120.0, 360.0 * protein_hue(2, 4));
    close(21.0 + 180.0 + 120.0, 360.0 * protein_hue(3, 4));
}

#[test]
fn scale_contig() {
    // ScaleContigTest
    let v = ScaleContig::new("1", &[200, 10]);
    assert_eq!(v.scale(10).unwrap(), 1);
    let first = 123_457;
    let v = ScaleContig::new("1", &[first + 1001, first, first + 1, first]);
    assert_eq!(v.scale(first).unwrap(), 1);
    assert_eq!(v.scale(first + 1).unwrap(), 1 + 10);
    assert_eq!(v.scale(first + 1001).unwrap(), 1 + 10 + 37);
    assert_eq!(log_distance(1), 10);
    assert_eq!(log_distance(10), 11);
    assert_eq!(log_distance(100), 18);
    assert_eq!(log_distance(1000), 37);
    let map = |m: &[(i32, i32)]| ScaleContig::from_map("Dummy", m.iter().copied().collect());
    let v = map(&[(1000, 10)]);
    assert_eq!((v.interpolate(900), v.interpolate(1000), v.interpolate(1100)), (10, 10, 10));
    let v = map(&[(1000, 10), (2000, 20)]);
    assert_eq!([900, 1000, 1100, 1500, 2000, 2100].map(|x| v.interpolate(x)), [10, 10, 11, 15, 20, 20]);
    let v = map(&[(1000, 10), (2000, 20), (3000, 40)]);
    assert_eq!([2000, 2100, 2500, 3000].map(|x| v.interpolate(x)), [20, 22, 30, 40]);
}

#[test]
fn scale_introns() {
    // ScaleIntronsTest
    let original = [100, 2000, 3001, 4600, 5601, 11000];
    let scaled = [100, 2000, 2011, 3610, 3621, 9020];
    let r = |s: i32, e: i32| Region::new("1", s, e);
    let exons = vec![r(original[0], original[1]), r(original[2], original[3]), r(original[4], original[5])];
    let v = ScaleIntrons::from_regions(&exons);
    assert_eq!(v.introns, vec![r(2001, 3000), r(4601, 5600)]);
    for i in 0..3 {
        assert_eq!(v.scale_region(&exons[i]), r(scaled[2 * i], scaled[2 * i + 1]));
    }
    assert_eq!(v.scale_region(&r(original[2] + 10, original[3] - 10)), r(scaled[2] + 10, scaled[3] - 10));
    assert_eq!(v.scale_region(&r(original[0] + 10, original[3] - 10)), r(scaled[0] + 10, scaled[3] - 10));
    assert_eq!(v.scale_region(&r(original[2] - 250, original[3] - 10)), r(scaled[2] - 3, scaled[3] - 10));
    assert_eq!(v.scale_region(&r(original[2] + 10, original[3] + 250)), r(scaled[2] + 10, scaled[3] + 3));
    assert_eq!(v.scale_region(&r(original[2] - 250, original[3] + 250)), r(scaled[2] - 3, scaled[3] + 3));
}

fn test_exon(chr: &str, start: i32, end: i32, rank: i32) -> VisGeneExon {
    VisGeneExon { cluster_id: 1, gene: chr.into(), transcript: chr.into(), chromosome: chr.into(), annotation_type: AnnotationType::EXON_LOST, exon_rank: rank, exon_start: start, exon_end: end }
}

fn test_fusion(strand_up: i32, pos_up: i32, exon_up: i32, strand_down: i32, pos_down: i32, exon_down: i32) -> VisFusion {
    VisFusion {
        cluster_id: 1,
        reportable: true,
        gene_name_up: "UP".into(),
        transcript_up: "UP".into(),
        chr_up: "UP".into(),
        pos_up,
        strand_up,
        region_type_up: "UP".into(),
        fused_exon_up: exon_up,
        gene_name_down: "DOWN".into(),
        transcript_down: "DOWN".into(),
        chr_down: "DOWN".into(),
        pos_down,
        strand_down,
        region_type_down: "DOWN".into(),
        fused_exon_down: exon_down,
    }
}

#[test]
fn disrupted_exons() {
    // DisruptedExonsTest
    let normal = |c: &str| vec![test_exon(c, 1, 100, 1), test_exon(c, 201, 300, 2)];
    let reverse = |c: &str| vec![test_exon(c, 1, 100, 2), test_exon(c, 201, 300, 1)];
    let regions = |exons: Vec<VisGeneExon>, f: VisFusion| -> Vec<(String, i32, i32)> {
        // disruptedGeneRegions(fusion, exons): the regions of the fusion only (no exon-lost entries here, as
        // the test calls the per-fusion method)
        let all = disrupted_gene_regions(&[f], &exons);
        let lost = exons.iter().filter(|e| e.annotation_type == AnnotationType::EXON_LOST).count();
        all[lost..].iter().map(|r| (r.chromosome.clone(), r.start, r.end)).collect()
    };
    let n: Vec<VisGeneExon> = normal("UP").into_iter().chain(normal("DOWN")).collect();
    assert_eq!(regions(n.clone(), test_fusion(1, 150, 1, 1, 150, 2)), vec![("UP".into(), 100, 300), ("DOWN".into(), 1, 201)]);
    assert_eq!(regions(n, test_fusion(1, 50, 1, 1, 250, 2)), vec![("UP".into(), 50, 300), ("DOWN".into(), 1, 250)]);
    let r: Vec<VisGeneExon> = reverse("UP").into_iter().chain(reverse("DOWN")).collect();
    assert_eq!(regions(r.clone(), test_fusion(-1, 150, 1, -1, 150, 2)), vec![("UP".into(), 1, 201), ("DOWN".into(), 100, 300)]);
    assert_eq!(regions(r, test_fusion(-1, 250, 1, -1, 50, 2)), vec![("UP".into(), 1, 250), ("DOWN".into(), 50, 300)]);
}

#[test]
fn fused_exons_convert_region() {
    // FusedExonsTest
    let gene = Region::new("1", 10000, 15000);
    assert_eq!(convert_region(1, &gene, 5000, 20000), (0, 5000));
}

#[test]
fn segments_tracks() {
    // SegmentsTest.testTrackStrategies
    let seg = |c: &str| {
        let mut s = VisSegment::new(1, 1, c, "1", "1000", 0.0, false).unwrap();
        s.track = 1;
        s
    };
    let r = increment_on_chromosome(vec![seg("1"), seg("2"), seg("1")], &[], false);
    assert_eq!(r.iter().map(|s| s.track).collect::<Vec<_>>(), vec![1, 1, 2]);
    assert_eq!(r.iter().map(|s| s.frame).collect::<Vec<_>>(), vec![1, 3, 5]);
}
