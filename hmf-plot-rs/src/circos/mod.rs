//! A Circos ring engine: draws what Circos 0.69 draws for a configuration file, as a display list.
//!
//! Supported (what the hmftools configurations use): karyotypes (hg19, hg38 built in, or a file); ideograms
//! with spacing, radius, thickness, fill, stroke and parallel labels; `<highlights>` (ideogram and radial);
//! plot tracks of type scatter (circle glyphs), histogram (joined bins, fill and outline), heatmap, highlight
//! and text (Circos's label snuggling, "array" method), each with `<axes>`, `<backgrounds>` and `<rules>`;
//! `<links>` as bezier curves or ribbons, with bezier radius, purity and crest.
//!
//! Geometry follows Circos.pm: positions are angles around a circle of the displayed ideograms and their
//! spacings; relative radii below 1 scale the ideogram's inner radius and from 1 up its outer radius; arcs
//! are polylines with a point every 0.5 degrees (dropping points closer than 10 pixels). What differs: shapes
//! are anti-aliased (GD fills polygons and draws thick lines without anti-aliasing), and text is rendered by
//! FreeType through BUI rather than by GD.

use crate::fmath::FMath;

pub mod conf;
pub mod expr;

use std::collections::HashMap;
use std::f64::consts::PI;

use crate::color::{CircosColours, Rgba};
use crate::display::{Figure, LineCap, Stroke};
use crate::png_out::PngMeta;
use crate::text::{self, Font};
use conf::{Block, Unit};

const DEG2RAD: f64 = PI / 180.0;
const RAD2DEG: f64 = 180.0 / PI;
/// housekeeping.conf
const ANGLE_STEP: f64 = 0.5;
const MIN_SLICE_STEP: f64 = 10.0;
const BEZIER_SAMPLES: usize = 40;
/// GD renders TrueType sizes in points at 96 dpi
const GD_DPI: f64 = 96.0;

pub const KARYOTYPE_HG38: &str = include_str!("../../resources/circos/karyotype.human.hg38.txt");
pub const KARYOTYPE_HG19: &str = include_str!("../../resources/circos/karyotype.human.hg19.txt");

/// Options that change what is drawn, for behaviours of a configuration that a port may correct.
#[derive(Clone, Debug, Default)]
pub struct Options {
    /// Ignore `skip_run` (Circos drops a data point whose value text equals the previous point's).
    pub ignore_skip_run: bool,
}

/// One displayed chromosome.
#[derive(Clone, Debug)]
struct Ideogram {
    chr: String,
    label: String,
    color: String,
    start: i64,
    end: i64,
    /// relative position (bases, spacings included) of its start around the circle
    rel_start: f64,
}

/// A data line: coordinates, value and options.
#[derive(Clone, Debug)]
struct Datum {
    chr: String,
    start: i64,
    end: i64,
    value: String,
    params: Vec<(String, String)>,
    hidden: bool,
}

impl Datum {
    fn param(&self, key: &str) -> Option<&str> {
        self.params.iter().rev().find(|(k, _)| k == key).map(|(_, v)| v.as_str())
    }

    fn num(&self) -> Option<f64> {
        self.value.trim().parse().ok()
    }
}

#[derive(Clone, Debug)]
struct Link {
    a: (String, i64, i64),
    b: (String, i64, i64),
    params: Vec<(String, String)>,
}

/// The defaults of Circos's etc/tracks/<type>.conf.
fn track_default(kind: &str, key: &str) -> Option<&'static str> {
    let table: &[(&str, &str)] = match kind {
        "scatter" => &[("glyph", "circle"), ("glyph_size", "10"), ("color", "grey"), ("stroke_color", "black"), ("stroke_thickness", "0"), ("r1", "0.79r"), ("r0", "0.7r"), ("orientation", "out")],
        "histogram" => &[("color", "black"), ("stroke_thickness", "1"), ("r1", "0.89r"), ("r0", "0.8r"), ("orientation", "out")],
        "heatmap" => &[("color", "spectral-11-div"), ("color_mapping", "0"), ("scale_log_base", "1"), ("stroke_thickness", "1"), ("r1", "0.975r"), ("r0", "0.95r")],
        "highlight" => &[("fill_color", "grey"), ("color", "black"), ("stroke_thickness", "1"), ("r1", "0.995r"), ("r0", "0.98r")],
        "text" => &[
            ("label_font", "default"),
            ("label_size", "12"),
            ("color", "black"),
            ("r0", "0.85r"),
            ("r1", "0.95r"),
            ("show_links", "no"),
            ("padding", "0p"),
            ("rpadding", "0p"),
            ("label_snuggle", "no"),
            ("max_snuggle_distance", "1r"),
            ("snuggle_sampling", "1"),
            ("snuggle_tolerance", "0.25r"),
            ("snuggle_refine", "no"),
        ],
        "link" => &[("ribbon", "no"), ("color", "black"), ("thickness", "1"), ("radius", "0.40r"), ("bezier_radius", "0r"), ("crest", "0.5"), ("bezier_radius_purity", "0.75")],
        "highlight.bg" => &[("fill_color", "grey_a4"), ("stroke_color", "black_a4"), ("stroke_thickness", "0"), ("r1", "1.2r"), ("r0", "0.2r"), ("ideogram", "0")],
        "axis" => &[("thickness", "1"), ("color", "grey")],
        _ => &[],
    };
    table.iter().find(|(k, _)| *k == key).map(|(_, v)| *v)
}

/// A parameter of a datum, else of its track, else the track type's default.
fn seek<'a>(datum: Option<&'a [(String, String)]>, track: &'a Block, kind: &str, key: &str) -> Option<&'a str> {
    if let Some(p) = datum {
        if let Some((_, v)) = p.iter().rev().find(|(k, _)| k == key) {
            return Some(v.as_str());
        }
    }
    track.get(key).or_else(|| track_default(kind, key))
}

fn perl_round(x: f64) -> f64 {
    expr::perl_round(x)
}

/// Reads data files by name (relative to the configuration's folder).
pub type Loader<'a> = dyn Fn(&str) -> Result<String, String> + 'a;

pub struct Circos<'a> {
    conf: Block,
    load: &'a Loader<'a>,
    opts: Options,
    colours: CircosColours,
    image_radius: f64,
    angle_offset: f64,
    ideograms: Vec<Ideogram>,
    /// the displayed chromosomes by name
    by_chr: HashMap<String, usize>,
    gcircum: f64,
    radius_outer: f64,
    radius_inner: f64,
    fig: Figure,
}

impl<'a> Circos<'a> {
    /// Lays out the ideograms of a configuration; `load` reads the files it names (karyotypes other than the
    /// built-in human ones, data files).
    pub fn new(conf_text: &str, load: &'a Loader<'a>, opts: Options) -> Result<Circos<'a>, String> {
        let conf = conf::parse(conf_text)?;
        // etc/image.conf, etc/image.generic.conf, etc/background.white.conf
        let image = conf.block("image");
        let image_radius = image.and_then(|b| b.get("radius")).and_then(Unit::parse).map(|u| u.value()).unwrap_or(1500.0);
        let angle_offset = image.and_then(|b| b.get("angle_offset")).and_then(|v| v.parse().ok()).unwrap_or(-90.0);
        let steps = image.and_then(|b| b.get("auto_alpha_steps")).and_then(|v| v.parse().ok()).unwrap_or(5);
        let mut colours = CircosColours::new(steps);
        for block in conf.blocks("colors") {
            for (k, v) in &block.params {
                colours.define(k, v);
            }
        }
        let background = image
            .and_then(|b| b.get("background"))
            .and_then(|c| colours.resolve(c))
            .unwrap_or(Rgba::WHITE);
        let size = (2.0 * image_radius).round() as u32;
        let mut fig = Figure::new(size, size, background);
        fig.meta = PngMeta { dpi: Some(GD_DPI), transparent_key: Some([0, 0, 1]) };

        let mut c = Circos {
            conf,
            load,
            opts,
            colours,
            image_radius,
            angle_offset,
            ideograms: Vec::new(),
            by_chr: HashMap::new(),
            gcircum: 0.0,
            radius_outer: 0.0,
            radius_inner: 0.0,
            fig,
        };
        c.layout()?;
        Ok(c)
    }

    fn layout(&mut self) -> Result<(), String> {
        let karyotype = self.conf.get("karyotype").ok_or("no karyotype")?.to_string();
        let text = if karyotype.ends_with("karyotype.human.hg38.txt") {
            KARYOTYPE_HG38.to_string()
        } else if karyotype.ends_with("karyotype.human.hg19.txt") {
            KARYOTYPE_HG19.to_string()
        } else {
            (self.load)(&karyotype)?
        };
        // chromosomes: "-hsY" hides, a list without '-' shows only those (with chromosomes_display_default no)
        let display_default = self.conf.flag("chromosomes_display_default").unwrap_or(true);
        let mut hidden = Vec::new();
        let mut shown = Vec::new();
        if let Some(list) = self.conf.get("chromosomes") {
            for item in list.split([';', ',']).map(|s| s.trim()).filter(|s| !s.is_empty()) {
                if let Some(h) = item.strip_prefix('-') {
                    hidden.push(h.to_string());
                } else {
                    shown.push(item.to_string());
                }
            }
        }
        for line in text.lines() {
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() < 7 || f[0] != "chr" {
                continue;
            }
            let chr = f[2].to_string();
            let display = (display_default || shown.contains(&chr)) && !hidden.contains(&chr);
            if !display {
                continue;
            }
            let start: i64 = f[4].parse().map_err(|_| format!("karyotype: bad start in '{line}'"))?;
            let end: i64 = f[5].parse().map_err(|_| format!("karyotype: bad end in '{line}'"))?;
            self.by_chr.insert(chr.clone(), self.ideograms.len());
            self.ideograms.push(Ideogram { chr, label: f[3].to_string(), color: f[6].to_string(), start, end, rel_start: 0.0 });
        }
        if self.ideograms.is_empty() {
            return Err("no ideograms to draw".into());
        }
        // sizes: each ideogram is a cover of scale 1 whose length is end - start; spacing is
        // "<x>r" of the displayed genome size (GSIZE_NOSCALE) or "<x>u" chromosome units
        let units: f64 = self.conf.get("chromosomes_units").and_then(|v| v.parse().ok()).unwrap_or(1.0);
        let gsize: f64 = self.ideograms.iter().map(|i| (i.end - i.start) as f64).sum();
        let ideo = self.conf.block("ideogram").cloned().unwrap_or_default();
        let spacing = match ideo.block("spacing").and_then(|s| s.get("default")).and_then(Unit::parse) {
            Some(Unit::R(v)) => v * gsize,
            Some(Unit::U(v)) => v * units,
            Some(u) => u.value(),
            None => 0.0,
        };
        let mut pos = 0.0;
        for i in &mut self.ideograms {
            i.rel_start = pos;
            pos += (i.end - i.start) as f64 + spacing;
        }
        self.gcircum = pos;
        let radius = match ideo.get("radius").and_then(Unit::parse) {
            Some(Unit::R(v)) => v * self.image_radius,
            Some(u) => u.value(),
            None => 0.9 * self.image_radius,
        };
        let thickness = match ideo.get("thickness").and_then(Unit::parse) {
            Some(Unit::R(v)) => v * self.image_radius,
            Some(u) => u.value(),
            None => 0.0,
        };
        self.radius_outer = radius;
        self.radius_inner = radius - thickness;
        Ok(())
    }

    // -----------------------------------------------------------------------------------------------------
    // geometry

    /// The ideogram holding a position (Circos's get_ideogram_idx).
    fn ideogram_of(&self, chr: &str, pos: i64) -> Option<usize> {
        let &i = self.by_chr.get(chr)?;
        let id = &self.ideograms[i];
        (pos >= id.start && pos <= id.end).then_some(i)
    }

    /// Circos's getanglepos: degrees, clockwise from 3 o'clock, angle_offset added.
    fn angle(&self, chr: &str, pos: f64) -> f64 {
        let rel = match self.by_chr.get(chr) {
            Some(&i) => self.ideograms[i].rel_start + (pos - self.ideograms[i].start as f64),
            None => 0.0,
        };
        360.0 / self.gcircum * rel + self.angle_offset
    }

    fn xy(&self, angle: f64, r: f64) -> (f64, f64) {
        (self.image_radius + r * (angle * DEG2RAD).fcos(), self.image_radius + r * (angle * DEG2RAD).fsin())
    }

    /// Circos's unit_parse of a radius: "<x>r" is x times the ideogram's inner radius when x < 1 and its
    /// outer radius otherwise; "<x>p" pixels; an expression with dims()/conf() is evaluated.
    fn radius(&self, v: &str) -> Result<f64, String> {
        match Unit::parse(v) {
            Some(Unit::R(x)) => Ok(x * if x < 1.0 { self.radius_inner } else { self.radius_outer }),
            Some(u) => Ok(u.value()),
            None => expr::eval(v, &DimsResolver { c: self }),
        }
    }

    fn colour(&self, spec: Option<&str>) -> Option<Rgba> {
        spec.and_then(|s| self.colours.resolve(s))
    }

    /// Circos's slice(): the polygon (or, for one radius, the arc) between two positions and two radii.
    fn slice_points(&self, chr: &str, start: f64, end: f64, r_from: f64, r_to: f64) -> (Vec<(f64, f64)>, bool) {
        let (mut a0, mut a1) = (self.angle(chr, start), self.angle(chr, end));
        if a1 < a0 {
            std::mem::swap(&mut a0, &mut a1);
        }
        if a0 > a1 {
            a0 -= 360.0;
        }
        let mut pts = Vec::new();
        let (mut xp, mut yp) = (0.0, 0.0);
        let mut a = a0;
        while a <= a1 {
            let (x, y) = self.xy(a, r_from);
            let d = ((x - xp).fpowi(2) + (y - yp).fpowi(2)).sqrt();
            a += ANGLE_STEP;
            if xp != 0.0 && yp != 0.0 && d < MIN_SLICE_STEP {
                continue;
            }
            pts.push((x, y));
            (xp, yp) = (x, y);
        }
        if a1 != a0 {
            pts.push(self.xy(a1, r_from));
        }
        let closed = r_from != r_to;
        if closed {
            let (mut xp, mut yp) = (0.0, 0.0);
            let mut a = a1;
            while a > a0 {
                let (x, y) = self.xy(a, r_to);
                let d = ((x - xp).fpowi(2) + (y - yp).fpowi(2)).sqrt();
                a -= ANGLE_STEP;
                if xp != 0.0 && yp != 0.0 && d < MIN_SLICE_STEP {
                    continue;
                }
                pts.push((x, y));
                (xp, yp) = (x, y);
            }
            pts.push(self.xy(a0, r_to));
        }
        (pts, closed)
    }

    /// Draws a slice: filled when `fill` is given, outlined when `stroke` > 0.
    fn slice(&mut self, chr: &str, start: f64, end: f64, r_from: f64, r_to: f64, fill: Option<Rgba>, edge: Option<Rgba>, edge_width: f64) {
        let (pts, closed) = self.slice_points(chr, start, end, r_from, r_to);
        let stroke = match edge {
            Some(c) if edge_width > 0.0 => Some(Stroke { color: c, width: edge_width, dash: Vec::new(), cap: LineCap::Round }),
            _ => None,
        };
        // a slice of no width is a radial line: GD fills the single pixel column of its polygon
        let degenerate = start == end;
        if degenerate {
            let (p, q) = (self.xy(self.angle(chr, start), r_from), self.xy(self.angle(chr, start), r_to));
            let s = stroke.or(fill.map(|c| Stroke { color: c, width: 1.0, dash: Vec::new(), cap: LineCap::Round }));
            if let Some(s) = s {
                self.fig.line(p.0, p.1, q.0, q.1, s);
            }
            return;
        }
        if closed {
            if fill.is_some() || stroke.is_some() {
                self.fig.polygon(pts, fill, stroke);
            }
        } else if let Some(s) = stroke {
            self.fig.polyline(pts, s);
        }
    }

    // -----------------------------------------------------------------------------------------------------
    // data

    fn read_data(&self, file: &str, link: bool, skip_run: bool) -> Result<(Vec<Datum>, Vec<Link>), String> {
        let text = (self.load)(file)?;
        let mut data = Vec::new();
        let mut links = Vec::new();
        let mut prev_value: Option<String> = None;
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split_whitespace().collect();
            if link {
                if f.len() < 6 {
                    continue;
                }
                if !self.by_chr.contains_key(f[0]) || !self.by_chr.contains_key(f[3]) {
                    continue;
                }
                let n = |s: &str| s.parse::<i64>().map_err(|_| format!("{file}: bad position in '{line}'"));
                links.push(Link {
                    a: (f[0].to_string(), n(f[1])?, n(f[2])?),
                    b: (f[3].to_string(), n(f[4])?, n(f[5])?),
                    params: f.get(6).map(|o| conf::parse_options(o)).unwrap_or_default(),
                });
                continue;
            }
            if f.len() < 3 || !self.by_chr.contains_key(f[0]) {
                continue;
            }
            let n = |s: &str| s.parse::<i64>().map_err(|_| format!("{file}: bad position in '{line}'"));
            let value = f.get(3).map(|s| s.to_string()).unwrap_or_default();
            if skip_run && !self.opts.ignore_skip_run && f.len() > 3 {
                if prev_value.as_deref() == Some(value.as_str()) {
                    continue;
                }
            }
            if f.len() > 3 {
                prev_value = Some(value.clone());
            }
            data.push(Datum {
                chr: f[0].to_string(),
                start: n(f[1])?,
                end: n(f[2])?,
                value,
                params: f.get(4).map(|o| conf::parse_options(o)).unwrap_or_default(),
                hidden: false,
            });
        }
        Ok((data, links))
    }

    /// Applies a track's `<rules>` to its data (Circos::Rule: the first rule that passes applies its
    /// parameters and stops, unless its flow is "continue"; `show = no` hides the point).
    fn apply_rules(&self, track: &Block, data: &mut [Datum]) -> Result<(), String> {
        let Some(rules) = track.block("rules") else { return Ok(()) };
        let rules: Vec<&Block> = rules.blocks("rule").collect();
        for d in data.iter_mut() {
            for rule in &rules {
                let Some(cond) = rule.get("condition") else { continue };
                let pass = expr::eval(cond, &DatumResolver { d }).map_err(|e| format!("rule '{cond}': {e}"))? != 0.0;
                if pass {
                    for (k, v) in &rule.params {
                        if k == "condition" || k == "flow" {
                            continue;
                        }
                        if k == "show" {
                            d.hidden = !conf::is_true(v);
                            continue;
                        }
                        let value = if v.starts_with("eval(") {
                            let x = expr::eval(v, &DatumResolver { d })?;
                            format_number(x)
                        } else {
                            v.clone()
                        };
                        d.params.push((k.clone(), value));
                    }
                }
                let flow_continue = rule.get("flow").is_some_and(|f| f.starts_with("continue"));
                if pass && !flow_continue {
                    break;
                }
            }
        }
        Ok(())
    }

    // -----------------------------------------------------------------------------------------------------
    // drawing

    /// Draws everything and returns the figure.
    pub fn draw(mut self) -> Result<Figure, String> {
        let highlights = self.read_highlights()?;
        let ideo_conf = self.conf.block("ideogram").cloned().unwrap_or_default();
        for i in 0..self.ideograms.len() {
            self.draw_highlights(&highlights, i, false);
            let id = self.ideograms[i].clone();
            let fill = if ideo_conf.flag("fill").unwrap_or(false) {
                self.colours.resolve(ideo_conf.get("fill_color").unwrap_or(&id.color)).or_else(|| self.colours.resolve(&id.color))
            } else {
                None
            };
            let (ri, ro) = (self.radius_inner, self.radius_outer);
            self.slice(&id.chr, id.start as f64, id.end as f64, ri, ro, fill, None, 0.0);
            if ideo_conf.flag("show_label").unwrap_or(false) {
                self.draw_ideogram_label(&ideo_conf, &id)?;
            }
            self.draw_highlights(&highlights, i, true);
            let edge = self.colour(ideo_conf.get("stroke_color"));
            let width = ideo_conf.get("stroke_thickness").and_then(Unit::parse).map(|u| u.value()).unwrap_or(0.0);
            if width > 0.0 {
                self.slice(&id.chr, id.start as f64, id.end as f64, ri, ro, None, edge, width);
            }
        }
        self.draw_links()?;
        self.draw_plots()?;
        Ok(self.fig)
    }

    fn read_highlights(&self) -> Result<Vec<(Block, Vec<Datum>)>, String> {
        let mut out = Vec::new();
        if let Some(hs) = self.conf.block("highlights") {
            for h in hs.blocks("highlight") {
                let Some(file) = h.get("file") else { continue };
                let (data, _) = self.read_data(file, false, false)?;
                out.push((h.clone(), data));
            }
        }
        // Circos sorts highlight tracks by z (stable)
        out.sort_by(|a, b| zval(&a.0).total_cmp(&zval(&b.0)));
        Ok(out)
    }

    /// Circos::Track::Highlight::draw_highlights for one ideogram: the tracks whose `ideogram` setting is
    /// `on_ideogram` (drawn after the ideogram's label) or not (drawn before the ideogram).
    fn draw_highlights(&mut self, tracks: &[(Block, Vec<Datum>)], ideo: usize, on_ideogram: bool) {
        let id = self.ideograms[ideo].clone();
        for (track, data) in tracks {
            for d in data {
                if d.chr != id.chr {
                    continue;
                }
                let (s, e) = (d.start.max(id.start), d.end.min(id.end));
                if s > e {
                    continue;
                }
                let p = Some(d.params.as_slice());
                let is_ideo = seek(p, track, "highlight.bg", "ideogram").map(conf::is_true).unwrap_or(false);
                if is_ideo != on_ideogram {
                    continue;
                }
                let (r0, r1) = if is_ideo {
                    (self.radius_inner, self.radius_outer)
                } else {
                    let r0 = seek(p, track, "highlight.bg", "r0").map(|v| self.radius(v)).transpose().ok().flatten().unwrap_or(0.0);
                    let r1 = seek(p, track, "highlight.bg", "r1").map(|v| self.radius(v)).transpose().ok().flatten().unwrap_or(0.0);
                    (r0, r1)
                };
                let fill = self.colour(seek(p, track, "highlight.bg", "fill_color"));
                let edge = self.colour(seek(p, track, "highlight.bg", "stroke_color"));
                let width = seek(p, track, "highlight.bg", "stroke_thickness").and_then(Unit::parse).map(|u| u.value()).unwrap_or(0.0);
                self.slice(&id.chr, s as f64, e as f64, r0, r1, fill, edge, width);
            }
        }
    }

    /// The label size Circos takes from GD's string box: ink width and height at the GD pixel size.
    fn label_size(font: Font, size_pt: f64, text: &str) -> (f64, f64) {
        let e = text::extents(font, size_pt * GD_DPI / 72.0, text);
        (e.right - e.left, e.ascent + e.descent.max(0.0))
    }

    fn draw_ideogram_label(&mut self, ideo: &Block, id: &Ideogram) -> Result<(), String> {
        let font = Font::circos(ideo.get("label_font").unwrap_or("default")).unwrap_or(Font::CircosBold);
        let size = ideo.get("label_size").and_then(Unit::parse).map(|u| u.value()).unwrap_or(12.0);
        let mut label = id.label.clone();
        if ideo.get("label_case") == Some("upper") {
            label = label.to_uppercase();
        } else if ideo.get("label_case") == Some("lower") {
            label = label.to_lowercase();
        }
        let (w, h) = Self::label_size(font, size, &label);
        let parallel = ideo.flag("label_parallel").unwrap_or(false);
        let mut radius = match ideo.get("label_radius") {
            Some(v) => self.radius(v)?,
            None => self.radius_outer,
        };
        if ideo.flag("label_center").unwrap_or(false) {
            radius -= if parallel { h } else { w } / 2.0;
        }
        let angle = self.angle(&id.chr, (id.start + id.end) as f64 / 2.0);
        let color = self.colour(ideo.get("label_color")).unwrap_or(Rgba::BLACK);
        let rotated = ideo.flag("label_rotated").unwrap_or(true);
        self.draw_text(&label, font, size, color, angle, radius, rotated, parallel);
        Ok(())
    }

    /// Circos::Text::draw_text: places a label at an angle and radius and draws it as GD's stringFT does.
    #[allow(clippy::too_many_arguments)]
    fn draw_text(&mut self, label: &str, font: Font, size_pt: f64, color: Rgba, angle: f64, radius: f64, rotated: bool, parallel: bool) {
        let q = angle_quadrant(angle);
        let text_angle = text_angle_2(angle, rotated, parallel);
        let (w, h) = Self::label_size(font, size_pt, label);
        let mut angle_offset = 0.0;
        if !parallel {
            let offset = RAD2DEG * h / 2.0 / radius;
            angle_offset += if rotated && q >= 2 { -offset } else { offset };
        } else {
            let offset = RAD2DEG * w / 2.0 / radius;
            angle_offset += if rotated && (q == 1 || q == 2) { offset } else { -offset };
        }
        let radius_offset = if !parallel {
            if rotated && q >= 2 { w } else { 0.0 }
        } else if rotated && (q == 1 || q == 2) {
            h
        } else {
            0.0
        };
        let (x, y) = self.xy(angle + angle_offset, radius + radius_offset);
        self.fig.text(x, y, label, font, size_pt * GD_DPI / 72.0, color, text_angle);
    }

    fn draw_links(&mut self) -> Result<(), String> {
        let Some(links_block) = self.conf.block("links").cloned() else { return Ok(()) };
        let mut tracks: Vec<Block> = links_block.blocks("link").cloned().collect();
        tracks.sort_by(|a, b| zval(a).total_cmp(&zval(b)));
        for track in &tracks {
            let Some(file) = track.get("file") else { continue };
            let (_, links) = self.read_data(file, true, false)?;
            for l in &links {
                let p = Some(l.params.as_slice());
                let (Some(i1), Some(i2)) = (self.ideogram_of(&l.a.0, l.a.1), self.ideogram_of(&l.b.0, l.b.1)) else { continue };
                let r1 = self.radius(seek(p, track, "link", "radius1").or(seek(p, track, "link", "radius")).unwrap_or("0.4r"))?;
                let r2 = self.radius(seek(p, track, "link", "radius2").or(seek(p, track, "link", "radius")).unwrap_or("0.4r"))?;
                let color = self.colour(seek(p, track, "link", "color")).unwrap_or(Rgba::BLACK);
                let bez = BezierParams {
                    bezier_radius: self.radius(seek(p, track, "link", "bezier_radius").unwrap_or("0r"))?,
                    purity: seek(p, track, "link", "bezier_radius_purity").and_then(|v| v.parse().ok()),
                    crest: seek(p, track, "link", "crest").and_then(|v| v.parse().ok()),
                };
                let (s1, e1) = (l.a.1.max(self.ideograms[i1].start), l.a.2.min(self.ideograms[i1].end));
                let (s2, e2) = (l.b.1.max(self.ideograms[i2].start), l.b.2.min(self.ideograms[i2].end));
                if seek(p, track, "link", "ribbon").map(conf::is_true).unwrap_or(false) {
                    // ribbon(): arc 1, bezier end1 -> end2, arc 2 backwards, bezier start2 -> start1
                    let mut pts = self.arc_points(&l.a.0, s1 as f64, e1 as f64, r1);
                    pts.extend(self.bezier(&l.a.0, e1 as f64, r1, &l.b.0, e2 as f64, r2, &bez));
                    pts.extend(self.arc_points(&l.b.0, e2 as f64, s2 as f64, r2));
                    pts.extend(self.bezier(&l.b.0, s2 as f64, r2, &l.a.0, s1 as f64, r1, &bez));
                    let edge = self.colour(seek(p, track, "link", "stroke_color"));
                    let edge_w = seek(p, track, "link", "stroke_thickness").and_then(Unit::parse).map(|u| u.value()).unwrap_or(0.0);
                    if s1 == e1 && s2 == e2 {
                        // a ribbon of single bases has no area: GD's anti-aliased fill of the doubled curve
                        // leaves a trace about 1.5 pixels wide
                        let curve = self.bezier(&l.a.0, e1 as f64, r1, &l.b.0, e2 as f64, r2, &bez);
                        self.fig.polyline(curve, Stroke { color, width: 1.5, dash: Vec::new(), cap: LineCap::Round });
                    } else {
                        let stroke = edge.filter(|_| edge_w > 0.0).map(|c| Stroke::solid(c, edge_w));
                        self.fig.polygon(pts, Some(color), stroke);
                    }
                } else {
                    let thickness = seek(p, track, "link", "thickness").and_then(Unit::parse).map(|u| u.value()).unwrap_or(1.0).round().max(1.0);
                    let m1 = (s1 + e1) as f64 / 2.0;
                    let m2 = (s2 + e2) as f64 / 2.0;
                    let curve = self.bezier(&l.a.0, m1, r1, &l.b.0, m2, r2, &bez);
                    self.fig.polyline(curve, Stroke { color, width: thickness, dash: Vec::new(), cap: LineCap::Round });
                }
            }
        }
        Ok(())
    }

    /// Circos's arc_points.
    fn arc_points(&self, chr: &str, start: f64, end: f64, radius: f64) -> Vec<(f64, f64)> {
        let (a0, a1) = (self.angle(chr, start), self.angle(chr, end));
        let step = if a0 < a1 { ANGLE_STEP } else { -ANGLE_STEP };
        let mut angles = Vec::new();
        let mut a = a0;
        if a0 < a1 {
            while a <= a1 {
                angles.push(a);
                a += step;
            }
        } else {
            while a >= a1 {
                angles.push(a);
                a += step;
                if a0 == a1 {
                    break;
                }
            }
        }
        let mut pts: Vec<(f64, f64)> = Vec::new();
        for a in angles {
            let (x, y) = self.xy(a, radius);
            if let Some(&(xp, yp)) = pts.last() {
                if ((x - xp).fpowi(2) + (y - yp).fpowi(2)).sqrt() < MIN_SLICE_STEP {
                    continue;
                }
            }
            pts.push((x, y));
            if a0 == a1 {
                break;
            }
        }
        pts.push(self.xy(a1, radius));
        pts
    }

    /// Circos's bezier_control_points and bezier_points: a curve between two positions through a control
    /// point at `bezier_radius` on the bisecting angle (moved by the purity), with crest points.
    #[allow(clippy::too_many_arguments)]
    fn bezier(&self, chr1: &str, pos1: f64, r1: f64, chr2: &str, pos2: f64, r2: f64, p: &BezierParams) -> Vec<(f64, f64)> {
        let (a1, a2) = (self.angle(chr1, pos1), self.angle(chr2, pos2));
        let (x1, y1) = self.xy(a1, r1);
        let (x2, y2) = self.xy(a2, r2);
        let size = 2.0 * self.image_radius;
        let bisecting = (((x1 + x2) / 2.0 - size / 2.0).fpowi(2) + ((y1 + y2) / 2.0 - size / 2.0).fpowi(2)).sqrt();
        let middle = if (a2 - a1).abs() > 180.0 { (a1 + a2 + 360.0) / 2.0 - 360.0 } else { (a2 + a1) / 2.0 };
        let mut br = p.bezier_radius;
        if let Some(k) = p.purity {
            let x = (1.0 - k).abs() * (br - bisecting).abs();
            br = if br > bisecting { if k > 1.0 { br + x } else { br - x } } else if k > 1.0 { br - x } else { br + x };
        }
        let mut ctrl = vec![(x1, y1), self.xy(middle, br), (x2, y2)];
        if let Some(crest) = p.crest {
            let cr = |r: f64| if r > br { r - (r - br).abs() * crest } else { r + (r - br).abs() * crest };
            ctrl.insert(1, self.xy(a1, cr(r1)));
            ctrl.insert(3, self.xy(a2, cr(r2)));
        }
        (0..BEZIER_SAMPLES)
            .map(|i| bezier_point(&ctrl, i as f64 / (BEZIER_SAMPLES - 1) as f64))
            .collect()
    }

    fn draw_plots(&mut self) -> Result<(), String> {
        let Some(plots) = self.conf.block("plots").cloned() else { return Ok(()) };
        let mut tracks: Vec<Block> = plots.blocks("plot").cloned().collect();
        tracks.sort_by(|a, b| zval(a).total_cmp(&zval(b)));
        for track in &tracks {
            let kind = track.get("type").unwrap_or("").to_string();
            let mut data = match track.get("file") {
                Some(f) => {
                    let skip_run = track.flag("skip_run").unwrap_or(false);
                    self.read_data(f, false, skip_run)?.0
                }
                None => Vec::new(),
            };
            self.apply_rules(track, &mut data)?;
            data.retain(|d| !d.hidden);
            // the track's value range: min and max when given, else the data's (calculate_track_statistics)
            let values: Vec<f64> = data.iter().filter_map(|d| d.num()).collect();
            let plot_min = track.get("min").and_then(|v| v.parse::<f64>().ok()).or_else(|| values.iter().copied().reduce(f64::min));
            let plot_max = track.get("max").and_then(|v| v.parse::<f64>().ok()).or_else(|| values.iter().copied().reduce(f64::max));
            self.draw_track_frame(track, &kind, plot_min, plot_max)?;
            match kind.as_str() {
                "scatter" => self.draw_scatter(track, data, plot_min, plot_max)?,
                "histogram" => self.draw_histogram(track, data, plot_min, plot_max)?,
                "heatmap" => self.draw_heatmap(track, data, plot_min, plot_max)?,
                "highlight" => self.draw_highlight_track(track, data)?,
                "text" => self.draw_text_track(track, data)?,
                _ => {}
            }
        }
        Ok(())
    }

    /// A track's `<backgrounds>` and `<axes>` on every ideogram.
    fn draw_track_frame(&mut self, track: &Block, kind: &str, plot_min: Option<f64>, plot_max: Option<f64>) -> Result<(), String> {
        let r0 = self.radius(seek(None, track, kind, "r0").unwrap_or("0r"))?;
        let r1 = self.radius(seek(None, track, kind, "r1").unwrap_or("0r"))?;
        let orientation_in = seek(None, track, kind, "orientation") == Some("in");
        for i in 0..self.ideograms.len() {
            let id = self.ideograms[i].clone();
            if let Some(bgs) = track.block("backgrounds") {
                for bg in bgs.blocks("background") {
                    let fill = self.colour(bg.get("color"));
                    let edge = self.colour(bg.get("stroke_color"));
                    let w = bg.get("stroke_thickness").and_then(Unit::parse).map(|u| u.value()).unwrap_or(0.0);
                    // y0 / y1 not used by the hmftools configurations: the whole track
                    self.slice(&id.chr, id.start as f64, id.end as f64, r0, r1, fill, edge, w);
                }
            }
            if let Some(axes) = track.block("axes") {
                let (min, max) = match (plot_min, plot_max) {
                    (Some(a), Some(b)) => (a, b),
                    _ => (r0, r1),
                };
                // Circos::Division::make_divisions: spacing ("<x>r" of the range), positions from min to max
                let mut divisions: Vec<(f64, f64, &Block)> = Vec::new();
                for axis in axes.blocks("axis") {
                    let spacing = match axis.get("spacing") {
                        Some(s) if s.ends_with('r') => s[..s.len() - 1].trim().parse::<f64>().unwrap_or(0.0) * (max - min),
                        Some(s) => s.parse::<f64>().unwrap_or(0.0),
                        None => 0.0,
                    };
                    if spacing <= 0.0 {
                        continue;
                    }
                    let tol = spacing / 1000.0;
                    let mut y = min;
                    while y <= max + tol {
                        divisions.push((spacing, y, axis));
                        y += spacing;
                    }
                }
                for (_, pos, axis) in divisions {
                    let color = self.colour(axis.get("color").or(track_default("axis", "color")));
                    let thickness = axis.get("thickness").or(track_default("axis", "thickness")).and_then(Unit::parse).map(|u| u.value()).unwrap_or(1.0);
                    let f = if max - min != 0.0 { (pos - min) / (max - min) } else { min };
                    let radius = if orientation_in { r1 - (r1 - r0) * f } else { r0 + (r1 - r0) * f };
                    self.slice(&id.chr, id.start as f64, id.end as f64, radius, radius, None, color, thickness);
                }
            }
        }
        Ok(())
    }

    /// Sorts data as Circos does for drawing: by z, then chromosome name (as text), then start; stable.
    fn sort_default(data: &mut [Datum]) {
        data.sort_by(|a, b| {
            let za = a.param("z").and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
            let zb = b.param("z").and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
            za.total_cmp(&zb).then_with(|| a.chr.cmp(&b.chr)).then_with(|| a.start.cmp(&b.start))
        });
    }

    /// check_data_range (trim to the ideogram) and get_ideogram_idx: the datum's ideogram and its trimmed
    /// coordinates, None when it lies outside every ideogram.
    fn place(&self, d: &Datum) -> Option<(usize, i64, i64)> {
        let &i = self.by_chr.get(&d.chr)?;
        let id = &self.ideograms[i];
        let (s, e) = (d.start.max(id.start), d.end.min(id.end));
        (s <= e).then_some((i, s, e))
    }

    /// check_value_limit: the value's radius and the radius of the value floor (0 when the range spans it).
    fn value_radius(value: f64, r0: f64, r1: f64, min: Option<f64>, max: Option<f64>, orientation_in: bool) -> (f64, f64) {
        let (Some(min), Some(max)) = (min, max) else { return (r1, r1) };
        let v = value.clamp(min, max);
        let floor = if min <= 0.0 && max >= 0.0 { 0.0 } else if min.abs() < max.abs() { min } else { max };
        let rd = (r1 - r0).abs();
        let dd = max - min;
        if dd == 0.0 {
            return (r1, r1);
        }
        if orientation_in {
            (r1 - rd * (v - min).abs() / dd, r1 - rd * (floor - min) / dd)
        } else {
            (r0 + rd * (v - min) / dd, r0 + rd * (floor - min) / dd)
        }
    }

    /// The colour of a datum from the track's colour list (encode_mapping, linear): one colour for all
    /// values, or a Brewer palette spread over [min, max).
    fn mapped_colour(&self, track_colour: &str, value: f64, min: Option<f64>, max: Option<f64>) -> Option<Rgba> {
        let list: Vec<String> = if self.colours.resolve(track_colour).is_some() {
            vec![track_colour.to_string()]
        } else if let Some(n) = track_colour.split('-').nth(1).and_then(|n| n.parse::<usize>().ok()) {
            (1..=n).map(|i| format!("{track_colour}-{i}")).collect()
        } else {
            track_colour.split(',').map(|s| s.trim().to_string()).collect()
        };
        if list.len() <= 1 {
            return self.colour(list.first().map(|s| s.as_str()));
        }
        let (min, max) = (min.unwrap_or(0.0), max.unwrap_or(1.0));
        let n = list.len();
        let i = if value < min {
            0
        } else if value >= max {
            n - 1
        } else {
            (((value - min) / (max - min)) * n as f64).floor() as usize
        };
        self.colour(Some(&list[i.min(n - 1)]))
    }

    fn draw_scatter(&mut self, track: &Block, mut data: Vec<Datum>, min: Option<f64>, max: Option<f64>) -> Result<(), String> {
        Self::sort_default(&mut data);
        let r0 = self.radius(seek(None, track, "scatter", "r0").unwrap_or("0r"))?;
        let r1 = self.radius(seek(None, track, "scatter", "r1").unwrap_or("0r"))?;
        let orientation_in = seek(None, track, "scatter", "orientation") == Some("in");
        let track_colour = seek(None, track, "scatter", "color").unwrap_or("grey").to_string();
        for d in &data {
            let Some((_, s, e)) = self.place(d) else { continue };
            let Some(value) = d.num() else { continue };
            let p = Some(d.params.as_slice());
            let size = seek(p, track, "scatter", "glyph_size").and_then(Unit::parse).map(|u| u.value()).unwrap_or(0.0);
            if size <= 0.0 {
                continue;
            }
            let colour = match d.param("fill_color").or(d.param("color")) {
                Some(c) => self.colour(Some(c)),
                None => self.mapped_colour(&track_colour, value, min, max),
            };
            let (radius, _) = Self::value_radius(value, r0, r1, min, max, orientation_in);
            let angle = self.angle(&d.chr, (s + e) as f64 / 2.0);
            let (x, y) = self.xy(angle, radius);
            let stroke_w = seek(p, track, "scatter", "stroke_thickness").and_then(Unit::parse).map(|u| u.value()).unwrap_or(0.0);
            let stroke = if stroke_w > 0.0 { self.colour(seek(p, track, "scatter", "stroke_color")).map(|c| Stroke::solid(c, stroke_w)) } else { None };
            // GD's filledArc of width w covers w + 1 pixels: radius (w + 1) / 2
            self.fig.circle(x, y, (size + 1.0) / 2.0, colour, stroke);
        }
        Ok(())
    }

    fn draw_histogram(&mut self, track: &Block, mut data: Vec<Datum>, min: Option<f64>, max: Option<f64>) -> Result<(), String> {
        Self::sort_default(&mut data);
        let r0 = self.radius(seek(None, track, "histogram", "r0").unwrap_or("0r"))?;
        let r1 = self.radius(seek(None, track, "histogram", "r1").unwrap_or("0r"))?;
        let orientation_in = seek(None, track, "histogram", "orientation") == Some("in");
        // (ideogram, start, end, radius, radius of the floor, fill, outline colour, outline width)
        let mut bins: Vec<(usize, f64, f64, f64, f64, Option<Rgba>, Option<Rgba>, f64)> = Vec::new();
        for d in &data {
            let Some((i, s, e)) = self.place(d) else { continue };
            let Some(value) = d.num() else { continue };
            let p = Some(d.params.as_slice());
            let (radius, radius0) = Self::value_radius(value, r0, r1, min, max, orientation_in);
            let fill = self.colour(seek(p, track, "histogram", "fill_color"));
            let colour = seek(p, track, "histogram", "color").and_then(|c| self.colours.resolve(c));
            let width = seek(p, track, "histogram", "thickness").or(seek(p, track, "histogram", "stroke_thickness")).and_then(Unit::parse).map(|u| u.value()).unwrap_or(0.0);
            bins.push((i, s as f64, e as f64, radius, radius0, fill, colour, width));
        }
        // A run of bins on one ideogram is joined (Circos's join_bins): each bin reaches to where the next
        // starts (the midpoint of a gap). Its fill is drawn as one polygon (the bins' outline over the floor),
        // which is the union of Circos's per-bin slices without anti-aliasing seams between them; its outline
        // is the bins' tops and the steps between them, closed down to the floor at both ends.
        let mut k = 0;
        while k < bins.len() {
            let i = bins[k].0;
            let mut end = k;
            while end + 1 < bins.len() && bins[end + 1].0 == i {
                end += 1;
            }
            let chr = self.ideograms[i].chr.clone();
            let run = &bins[k..=end];
            let mut starts: Vec<f64> = Vec::with_capacity(run.len());
            for (j, b) in run.iter().enumerate() {
                if j == 0 {
                    starts.push(b.1);
                } else {
                    let prev_end = run[j - 1].2;
                    starts.push(if prev_end == b.1 - 1.0 { b.1 } else { (prev_end + b.1) / 2.0 });
                }
            }
            let ends: Vec<f64> = (0..run.len()).map(|j| if j + 1 < run.len() { starts[j + 1] } else { run[j].2 }).collect();
            let radius0 = run[0].4;
            let mut top: Vec<(f64, f64)> = Vec::new();
            for (j, b) in run.iter().enumerate() {
                top.extend(self.arc_points(&chr, starts[j], ends[j], b.3));
            }
            if let Some(fill) = run[0].5 {
                let mut poly = top.clone();
                let mut floor = self.arc_points(&chr, starts[0], ends[run.len() - 1], radius0);
                floor.reverse();
                poly.extend(floor);
                self.fig.polygon(poly, Some(fill), None);
            }
            let (colour, width) = (run[0].6, run[0].7);
            if let Some(c) = colour.filter(|_| width > 0.0) {
                let mut outline = vec![self.xy(self.angle(&chr, starts[0]), radius0)];
                outline.extend(top);
                outline.push(self.xy(self.angle(&chr, ends[run.len() - 1]), radius0));
                self.fig.polyline(outline, Stroke { color: c, width, dash: Vec::new(), cap: LineCap::Round });
            }
            k = end + 1;
        }
        Ok(())
    }

    fn draw_heatmap(&mut self, track: &Block, mut data: Vec<Datum>, min: Option<f64>, max: Option<f64>) -> Result<(), String> {
        // Circos sorts heat maps by z, then by size, largest first
        data.sort_by(|a, b| (b.end - b.start).cmp(&(a.end - a.start)));
        let r0 = self.radius(seek(None, track, "heatmap", "r0").unwrap_or("0r"))?;
        let r1 = self.radius(seek(None, track, "heatmap", "r1").unwrap_or("0r"))?;
        let track_colour = seek(None, track, "heatmap", "color").unwrap_or("spectral-11-div").to_string();
        for d in &data {
            let Some((i, s, e)) = self.place(d) else { continue };
            let chr = self.ideograms[i].chr.clone();
            let value = d.num().unwrap_or(0.0);
            let colour = match d.param("color") {
                Some(c) => self.colour(Some(c)),
                None => self.mapped_colour(&track_colour, value, min, max),
            };
            let p = Some(d.params.as_slice());
            let edge = self.colour(d.param("stroke_color")).or(colour);
            let width = seek(p, track, "heatmap", "stroke_thickness").and_then(Unit::parse).map(|u| u.value()).unwrap_or(0.0);
            self.slice(&chr, s as f64, e as f64, r0, r1, colour, edge, width);
        }
        Ok(())
    }

    fn draw_highlight_track(&mut self, track: &Block, mut data: Vec<Datum>) -> Result<(), String> {
        Self::sort_default(&mut data);
        for d in &data {
            let Some((i, s, e)) = self.place(d) else { continue };
            let chr = self.ideograms[i].chr.clone();
            let p = Some(d.params.as_slice());
            let r0 = self.radius(seek(p, track, "highlight", "r0").unwrap_or("0r"))?;
            let r1 = self.radius(seek(p, track, "highlight", "r1").unwrap_or("0r"))?;
            let edge = self.colour(seek(p, track, "highlight", "stroke_color").or(seek(p, track, "highlight", "color")));
            let width = seek(p, track, "highlight", "stroke_thickness").and_then(Unit::parse).map(|u| u.value()).unwrap_or(0.0);
            let fill = self.colour(seek(p, track, "highlight", "fill_color"));
            self.slice(&chr, s as f64, e as f64, r0, r1, fill, edge, width);
        }
        Ok(())
    }

    /// A text track: Circos's label placement (layers, snuggling by pixel shifts, "array" stack method),
    /// then draw_text for every placed label.
    fn draw_text_track(&mut self, track: &Block, data: Vec<Datum>) -> Result<(), String> {
        let font = Font::circos(seek(None, track, "text", "label_font").unwrap_or("default")).unwrap_or(Font::CircosCondensed);
        let r1 = perl_round(self.radius(seek(None, track, "text", "r1").unwrap_or("0r"))?);
        let res = 1.0 * r1 * 2.0 * PI / 360.0; // text_pixel_subsampling = 1
        let snuggle = seek(None, track, "text", "label_snuggle").map(conf::is_true).unwrap_or(false);
        let mut labels: Vec<TextLabel> = Vec::new();
        for d in &data {
            let Some((_, s, e)) = self.place(d) else { continue };
            let p = Some(d.params.as_slice());
            let size = seek(p, track, "text", "label_size").and_then(Unit::parse).map(|u| u.value()).unwrap_or(12.0);
            let (w, h) = Self::label_size(font, size, &d.value);
            let radius = self.radius(seek(p, track, "text", "r0").unwrap_or("0r"))?;
            let rpadding = match seek(p, track, "text", "rpadding").and_then(Unit::parse) {
                Some(Unit::R(x)) => x * w,
                Some(u) => u.value(),
                None => 0.0,
            };
            let colour = self.colour(seek(p, track, "text", "color")).unwrap_or(Rgba::BLACK);
            labels.push(TextLabel {
                text: d.value.clone(),
                size,
                colour,
                angle: self.angle(&d.chr, (s + e) as f64 / 2.0),
                radius,
                dimr: w,
                dima: h,
                rpadding,
                layer: None,
                skip: false,
                angle_new: None,
                radius_shift: 0.0,
            });
        }
        // stable sort by angle
        labels.sort_by(|a, b| a.angle.total_cmp(&b.angle));
        let tolerance_r = match seek(None, track, "text", "snuggle_tolerance").and_then(Unit::parse) {
            Some(u) => u,
            None => Unit::N(0.0),
        };
        place_labels(&mut labels, res, r1, snuggle, self.angle_offset, tolerance_r);
        let rotated = track.flag("label_rotate").unwrap_or(true);
        let parallel = track.flag("label_parallel").unwrap_or(false);
        for l in &labels {
            if l.layer.is_none() {
                continue;
            }
            let angle = l.angle_new.unwrap_or(l.angle);
            let radius = l.radius + l.radius_shift + l.rpadding;
            self.draw_text(&l.text, font, l.size, l.colour, angle, radius, rotated, parallel);
        }
        Ok(())
    }
}

struct BezierParams {
    bezier_radius: f64,
    purity: Option<f64>,
    crest: Option<f64>,
}

/// A point of the Bezier curve of any degree (Math::Bezier).
fn bezier_point(ctrl: &[(f64, f64)], t: f64) -> (f64, f64) {
    let n = ctrl.len() - 1;
    let (mut x, mut y) = (0.0, 0.0);
    let mut binom = 1.0;
    for (k, &(cx, cy)) in ctrl.iter().enumerate() {
        if k > 0 {
            binom = binom * (n - k + 1) as f64 / k as f64;
        }
        let b = binom * t.fpowi(k as i32) * (1.0 - t).fpowi((n - k) as i32);
        x += b * cx;
        y += b * cy;
    }
    (x, y)
}

fn zval(b: &Block) -> f64 {
    b.get("z").and_then(|v| v.parse().ok()).unwrap_or(0.0)
}

fn format_number(x: f64) -> String {
    if x.fract() == 0.0 { format!("{}", x as i64) } else { format!("{x}") }
}

/// Circos::Geometry::angle_quadrant.
fn angle_quadrant(angle: f64) -> u8 {
    let mut a = angle;
    while a < -90.0 {
        a += 360.0;
    }
    while a > 270.0 {
        a -= 360.0;
    }
    if a <= 0.0 {
        0
    } else if a <= 90.0 {
        1
    } else if a <= 180.0 {
        2
    } else {
        3
    }
}

/// Circos::Text::text_angle_2: GD's text angle (degrees, counter-clockwise).
fn text_angle_2(angle: f64, rotated: bool, parallel: bool) -> f64 {
    let q = angle_quadrant(angle);
    let mut t = if q == 0 { -angle } else { 360.0 - angle };
    if rotated && q >= 2 {
        t += 180.0;
    }
    if parallel {
        if rotated {
            if q == 1 || q == 3 {
                t += 90.0;
            } else {
                t -= 90.0;
            }
        } else {
            t -= 90.0;
        }
    }
    t
}

#[derive(Clone, Debug)]
struct TextLabel {
    text: String,
    size: f64,
    colour: Rgba,
    angle: f64,
    radius: f64,
    /// radial extent (the label's width), angular extent (its height)
    dimr: f64,
    dima: f64,
    rpadding: f64,
    layer: Option<usize>,
    skip: bool,
    angle_new: Option<f64>,
    radius_shift: f64,
}

/// Circos's text track placement, "array" snuggle method, without refinement or link overlap tests (off by
/// default): passes alternate between finding the lowest height any label can take and placing the labels
/// that can take it, shifting each by up to `max_snuggle_distance` (1r: its height) pixels.
fn place_labels(labels: &mut [TextLabel], res: f64, r1: f64, snuggle: bool, angle_offset: f64, tolerance: Unit) {
    let n = (2.0 * 360.0 * res) as usize + 2;
    let mut stack = vec![0.0f64; n];
    let array_offset = 45.0;
    let mut layer = 0usize;
    let mut seek_min = true;
    let mut global_min: Option<f64> = Some(0.0);
    let mut label_not_placed = 0;
    let mut all_placed = 0usize;
    let mut iters: HashMap<usize, usize> = HashMap::new();
    let idx = |x: f64| -> usize { (x.trunc().max(0.0) as usize).min(n - 1) };
    loop {
        let mut label_placed = 0;
        for l in labels.iter_mut() {
            if l.layer.is_some() {
                continue;
            }
            if l.skip {
                l.skip = false;
                continue;
            }
            let range: Vec<f64> = if !snuggle {
                vec![0.0]
            } else {
                let sd_max = l.dima; // 1r of the label's height
                let ss = 1.0;
                let mut r = Vec::new();
                for k in 0..=((sd_max / ss) as i64) {
                    r.push(perl_round(-(k as f64) * ss));
                    r.push(perl_round(k as f64 * ss));
                }
                r.sort_by(|a, b| a.abs().total_cmp(&b.abs()));
                r
            };
            let mut label_min: Option<f64> = None;
            let mut best_shift = 0.0;
            let mult = RAD2DEG / l.radius;
            for &shift in &range {
                let angle_new = l.angle + mult * shift;
                let a = angle_new + array_offset - angle_offset;
                let h = stack[idx(a * res)];
                let ah = RAD2DEG * l.dima / (l.radius + h);
                let (a1, a2) = (a - ah / 2.0, a + ah / 2.0);
                let mut cur = 0.0f64;
                let (e0, e1) = (idx(a1 * res), idx(a2 * res));
                for e in e0..=e1.max(e0) {
                    cur = cur.max(stack[e]);
                }
                if l.radius + cur + l.dimr > r1 {
                    continue;
                }
                let d = cur - global_min.unwrap_or(0.0);
                let pass = if !seek_min {
                    let tol = match tolerance {
                        Unit::R(x) => x * l.dimr,
                        u => u.value(),
                    };
                    match label_min {
                        None => d <= tol,
                        Some(_) => d >= 0.0 && d <= tol && shift.abs() < f64::abs(best_shift),
                    }
                } else {
                    label_min.is_none_or(|m| cur < m)
                };
                if pass {
                    label_min = Some(cur);
                    if !seek_min {
                        l.angle_new = Some(angle_new);
                        best_shift = shift;
                    }
                }
            }
            if seek_min {
                let lm = label_min.unwrap_or(0.0);
                match global_min {
                    None => global_min = Some(lm),
                    Some(g) if lm - g < 0.0 => global_min = Some(lm),
                    Some(g) if lm - g > 0.0 => l.skip = true,
                    _ => {}
                }
                continue;
            }
            let Some(angle_new) = l.angle_new else { continue };
            let lm = label_min.unwrap_or(0.0);
            l.radius_shift = lm;
            l.layer = Some(layer);
            label_placed += 1;
            all_placed += 1;
            let ah_outer = RAD2DEG * l.dima / (l.radius + l.radius_shift + l.dimr);
            let span = |a: f64| perl_round(perl_round((a - angle_offset) * res));
            let (s0, s1) = (span(angle_new - ah_outer / 2.0), span(angle_new + ah_outer / 2.0));
            let (low, high) = (l.radius_shift + l.rpadding, l.radius_shift + l.rpadding + l.dimr);
            let mut a = s0 as i64;
            while a <= s1 as i64 {
                let i = idx(a as f64 + array_offset * res);
                stack[i] = stack[i].max(low).max(high);
                a += 1;
            }
        }
        if seek_min {
            seek_min = false;
        } else {
            seek_min = true;
            if label_placed == 0 {
                label_not_placed += 1;
                layer += 1;
                global_min = None;
            } else {
                label_not_placed = 0;
            }
            let c = iters.entry(all_placed).or_insert(0);
            *c += 1;
            if *c > 21 {
                break;
            }
        }
        if !(label_placed > 0 || label_not_placed < 2) {
            break;
        }
    }
}

struct DimsResolver<'b, 'a> {
    c: &'b Circos<'a>,
}

impl expr::Resolver for DimsResolver<'_, '_> {
    fn var(&self, _: &str) -> Option<f64> {
        None
    }
    fn dims(&self, args: &[&str]) -> Option<f64> {
        match args.last()? {
            &"radius_outer" | &"radius" => Some(self.c.radius_outer),
            &"radius_inner" => Some(self.c.radius_inner),
            &"radius_middle" => Some((self.c.radius_inner + self.c.radius_outer) / 2.0),
            &"thickness" => Some(self.c.radius_outer - self.c.radius_inner),
            _ => None,
        }
    }
    fn conf(&self, args: &[&str]) -> Option<f64> {
        let mut b = &self.c.conf;
        for a in &args[..args.len() - 1] {
            b = b.block(a)?;
        }
        Unit::parse(b.get(args.last()?)?).map(|u| u.value())
    }
}

struct DatumResolver<'a> {
    d: &'a Datum,
}

impl expr::Resolver for DatumResolver<'_> {
    fn var(&self, name: &str) -> Option<f64> {
        match name {
            "value" => self.d.num(),
            "start" => Some(self.d.start as f64),
            "end" => Some(self.d.end as f64),
            "size" => Some((self.d.end - self.d.start) as f64),
            _ => self.d.param(name).and_then(|v| v.parse().ok()),
        }
    }
    fn dims(&self, _: &[&str]) -> Option<f64> {
        None
    }
    fn conf(&self, _: &[&str]) -> Option<f64> {
        None
    }
}

/// Draws a Circos configuration file; data files are read relative to its folder (as Circos is run with
/// the configuration's folder as working directory by hmftools).
pub fn draw_conf_file(conf_path: &std::path::Path, opts: Options) -> Result<Figure, String> {
    let text = std::fs::read_to_string(conf_path).map_err(|e| format!("{}: {e}", conf_path.display()))?;
    let dir = conf_path.parent().map(|p| p.to_path_buf()).unwrap_or_default();
    let load = move |name: &str| -> Result<String, String> {
        let p = std::path::Path::new(name);
        let path = if p.is_absolute() { p.to_path_buf() } else { dir.join(p) };
        std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))
    };
    Circos::new(&text, &load, opts)?.draw()
}

#[cfg(test)]
mod tests;
