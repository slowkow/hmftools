//! Display-list tests of the Circos engine on small synthetic configurations.

use super::*;
use crate::display::Item;

const KARYO: &str = "chr - c1 1 0 100 chr1\nchr - c2 2 0 300 chr2\nchr - c3 3 0 50 chr3\n";

fn files(name: &str) -> Result<String, String> {
    match name {
        "k.txt" => Ok(KARYO.to_string()),
        "s.txt" => Ok("#chromosome\tstart\tend\tvalue\nc1\t50\t50\t0.5\nc2\t10\t10\t2.0\tcolor=255,0,0\nc2\t20\t20\t2.0\tcolor=255,0,0\n".into()),
        "h.txt" => Ok("c2\t0\t99\t1\nc2\t100\t199\t3\nc2\t200\t300\t-1\n".into()),
        "t.txt" => Ok("c2\t100\t100\tGENEA\tcolor=black\nc2\t100\t100\tGENEB\tcolor=black\n".into()),
        "l.txt" => Ok("c1\t50\t50\tc2\t150\t150\tcolor=blue\n".into()),
        _ => Err(format!("no file {name}")),
    }
}

fn conf(body: &str) -> String {
    format!(
        "karyotype = k.txt\nchromosomes = -c3\n<image>\nradius = 500p\nangle_offset = -90\n</image>\n<ideogram>\n<spacing>\ndefault = 0.1r\n</spacing>\nradius = 0.9r\nthickness = 50p\nfill = yes\nstroke_thickness = 0\n</ideogram>\n<colors>\nchr1 = 10,20,30\nchr2 = 40,50,60\n</colors>\n{body}"
    )
}

#[test]
fn layout_angles_and_radii() {
    let load = files;
    let c = Circos::new(&conf(""), &load, Options::default()).unwrap();
    // c3 hidden; spacing 0.1 of 400 bases = 40; circumference 100 + 40 + 300 + 40 = 480
    assert_eq!(c.ideograms.len(), 2);
    assert_eq!(c.gcircum, 480.0);
    assert_eq!(c.angle("c1", 0.0), -90.0);
    assert_eq!(c.angle("c2", 0.0), 360.0 * 140.0 / 480.0 - 90.0);
    assert_eq!((c.radius_outer, c.radius_inner), (450.0, 400.0));
    // relative radii below 1 scale the inner radius, from 1 up the outer one
    assert_eq!(c.radius("0.5r").unwrap(), 200.0);
    assert_eq!(c.radius("1.2r").unwrap(), 540.0);
    assert_eq!(c.radius("(dims(ideogram,radius_outer)+dims(ideogram,radius_inner))/2").unwrap(), 425.0);
    let fig = c.draw().unwrap();
    assert_eq!((fig.width, fig.height), (1000, 1000));
    // two ideogram polygons, filled with their karyotype colours
    let fills: Vec<String> = fig
        .items
        .iter()
        .filter_map(|i| match i {
            Item::Polygon { fill: Some(f), .. } => Some(f.to_hex()),
            _ => None,
        })
        .collect();
    assert_eq!(fills, ["#0A141E", "#28323C"]);
}

#[test]
fn scatter_radius_colour_and_skip_run() {
    let body = "<plots>\n<plot>\ntype = scatter\nfile = s.txt\nr0 = 0.5r\nr1 = 0.9r\nmin = 0\nmax = 1\ncolor = black\nglyph_size = 8\nskip_run = yes\n</plot>\n</plots>\n";
    let load = files;
    let circles = |opts: Options| -> Vec<(f64, f64, f64, String)> {
        let fig = Circos::new(&conf(body), &load, opts).unwrap().draw().unwrap();
        fig.items
            .iter()
            .filter_map(|i| match i {
                Item::Circle { cx, cy, r, fill: Some(f), .. } => Some((*cx, *cy, *r, f.to_hex())),
                _ => None,
            })
            .collect()
    };
    let with_skip = circles(Options::default());
    // the second "2.0" repeats the previous value: skip_run drops it
    assert_eq!(with_skip.len(), 2);
    // c1:50 is at -90 + 360 * 50 / 480 = -52.5 degrees; value 0.5 of [0, 1] between 200 and 360 px: 280
    let (cx, cy, r, ref colour) = with_skip[0];
    let a = (-52.5f64).to_radians();
    assert!((cx - (500.0 + 280.0 * a.cos())).abs() < 1e-9 && (cy - (500.0 + 280.0 * a.sin())).abs() < 1e-9);
    assert_eq!((r, colour.as_str()), (4.5, "#000000"));
    // values above max are drawn at max; the datum's colour wins over the track's
    assert_eq!(with_skip[1].3, "#FF0000");
    let all = circles(Options { ignore_skip_run: true });
    assert_eq!(all.len(), 3);
}

#[test]
fn rules_hide_and_eval() {
    let body = "<plots>\n<plot>\ntype = scatter\nfile = s.txt\nr0 = 0.5r\nr1 = 0.5r\nmin = 0\nmax = 5\nglyph_size = 8\ncolor = black\n<rules>\n<rule>\ncondition = var(value) < 1\nshow = no\nflow = continue\n</rule>\n<rule>\ncondition = var(value) > 1\nglyph_size = eval(int(remap_int(var(value),1,3,10,30)))\n</rule>\n</rules>\n</plot>\n</plots>\n";
    let load = files;
    let fig = Circos::new(&conf(body), &load, Options::default()).unwrap().draw().unwrap();
    let radii: Vec<f64> = fig.items.iter().filter_map(|i| if let Item::Circle { r, .. } = i { Some(*r) } else { None }).collect();
    // 0.5 hidden; the two 2.0 points get remap(2, 1, 3, 10, 30) = 20 -> radius 10.5 (GD: 21 pixels)
    assert_eq!(radii, [10.5, 10.5]);
}

#[test]
fn histogram_joins_bins_over_the_floor() {
    let body = "<plots>\n<plot>\ntype = histogram\nfile = h.txt\nr0 = 0.5r\nr1 = 0.75r\nmin = 0\nmax = 2\ncolor = black\nfill_color = red\n<axes>\n<axis>\ncolor = grey\nthickness = 1\nspacing = 0.5r\n</axis>\n</axes>\n<backgrounds>\n<background>\ncolor = vlgrey\n</background>\n</backgrounds>\n</plot>\n</plots>\n";
    let load = files;
    let fig = Circos::new(&conf(body), &load, Options::default()).unwrap().draw().unwrap();
    let red = CircosColours::new(5).resolve("red").unwrap();
    let fills: Vec<&Vec<(f64, f64)>> = fig
        .items
        .iter()
        .filter_map(|i| match i {
            Item::Polygon { points, fill: Some(f), .. } if *f == red => Some(points),
            _ => None,
        })
        .collect();
    // one polygon for the run of three bins on c2
    assert_eq!(fills.len(), 1);
    // the radii it reaches: floor 200 (value 0), 250 (1), 300 (3, clamped to 2), 200 (-1, clamped to 0)
    let radii: Vec<f64> = fills[0].iter().map(|&(x, y)| ((x - 500.0).powi(2) + (y - 500.0).powi(2)).sqrt()).collect();
    for r in &radii {
        assert!([200.0, 250.0, 300.0].iter().any(|v| (r - v).abs() < 1e-6), "radius {r}");
    }
    // axes at 0, 1 and 2 (spacing 0.5 of the range) on both ideograms: open arcs at 200, 250 and 300 px
    let axes = fig.items.iter().filter(|i| matches!(i, Item::Polyline { stroke, .. } if stroke.width == 1.0 && stroke.color == CircosColours::new(5).resolve("grey").unwrap())).count();
    assert_eq!(axes, 6);
}

#[test]
fn text_labels_snuggle_apart() {
    let body = "<plots>\n<plot>\ntype = text\nfile = t.txt\nr0 = 1.01r\nr1 = 1.5r\nlabel_size = 20p\nlabel_font = condensed\nlabel_snuggle = yes\n</plot>\n</plots>\n";
    let load = files;
    let fig = Circos::new(&conf(body), &load, Options::default()).unwrap().draw().unwrap();
    let texts: Vec<(f64, f64, f64)> = fig.items.iter().filter_map(|i| if let Item::Text { x, y, angle, .. } = i { Some((*x, *y, *angle)) } else { None }).collect();
    assert_eq!(texts.len(), 2);
    // the two labels share a position: the second is shifted by up to its height, or moved out a layer
    let d = ((texts[0].0 - texts[1].0).powi(2) + (texts[0].1 - texts[1].1).powi(2)).sqrt();
    assert!(d > 5.0, "labels overlap: {texts:?}");
}

#[test]
fn ribbon_of_single_bases_is_a_curve() {
    let body = "<links>\n<link>\nfile = l.txt\nradius = 0.5r\nbezier_radius = 0.1r\nribbon = yes\n</link>\n</links>\n";
    let load = files;
    let fig = Circos::new(&conf(body), &load, Options::default()).unwrap().draw().unwrap();
    let curves: Vec<&Vec<(f64, f64)>> = fig.items.iter().filter_map(|i| if let Item::Polyline { points, .. } = i { Some(points) } else { None }).collect();
    assert_eq!(curves.len(), 1);
    assert_eq!(curves[0].len(), BEZIER_SAMPLES);
    // it starts and ends on the link radius (200 px)
    let r = |p: (f64, f64)| ((p.0 - 500.0).powi(2) + (p.1 - 500.0).powi(2)).sqrt();
    assert!((r(curves[0][0]) - 200.0).abs() < 1e-9 && (r(*curves[0].last().unwrap()) - 200.0).abs() < 1e-9);
}

#[test]
fn quadrants_and_text_angles() {
    assert_eq!(angle_quadrant(-45.0), 0);
    assert_eq!(angle_quadrant(45.0), 1);
    assert_eq!(angle_quadrant(135.0), 2);
    assert_eq!(angle_quadrant(225.0), 3);
    // parallel, rotated labels read along the circle
    assert_eq!(text_angle_2(-75.0, true, true), 75.0 - 90.0);
    assert_eq!(text_angle_2(220.0, true, false), 360.0 - 220.0 + 180.0);
}
