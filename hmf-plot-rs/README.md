# hmf-plot-rs

Native drawing for the hmftools Rust ports (crate `hmf_plot`): the plots that the Java tools draw with R
(ggplot2) and Circos, drawn in Rust without R, Perl or Circos, and PDF documents. PURPLE's nine plots are the
first user of the raster side (`purple-rs/src/plots.rs`), LINX's visualiser the second
(`linx-rs/src/visualiser/plot.rs`, `fusion.rs`), and ORANGE's report the first user of the PDF side
(`orange-rs/src/report/`). CUPPA's summary plot (`cuppa-rs/src/plot.rs`) and qsee's report (`qsee-rs/src/plot/`)
use the ggplot2 model (`ggp`, `gtable`, `patchwork`) with both outputs: PNG through BUI and R-style PDF
(`pdf_figure`).

The raster renderer is [BUI](https://github.com/ppmpreetham/BUI) (Blender's UI kit as a C library), used only for
its software raster and its FreeType text engine. The PDF writer is pure Rust.

## Features

| Feature | What it gives | Needs |
|---|---|---|
| `raster` (default) | every module below but `pdf`: display lists, BUI's raster, PNG output, text, `gg`, `circos` | BUI's C kit: CMake and a C/C++ compiler |
| `pdf` | the `pdf` module (with `color`, `fmath` and `png_out`, which need no BUI) | nothing beyond Rust: pdf-writer, subsetter, ttf-parser, miniz_oxide |

A tool that only writes PDFs depends on `hmf-plot-rs = { path = "../hmf-plot-rs", default-features = false,
features = ["pdf"] }` (as orange-rs does) and builds without CMake. Default builds (purple-rs) are unchanged.

## Layers

| Module | What it is |
|---|---|
| `display` | A `Figure`: pixel size, background, PNG metadata, and a display list of `Item`s (rectangle, polygon, polyline with width, dash and cap, circle, text with rotation, clip) in device pixels. Nothing in it depends on BUI. |
| `raster` | Renders a display list with BUI's raster into RGB pixels. Adds what BUI lacks: butt caps and dashes, translucent multi-segment strokes drawn once (an offscreen BUI raster composited as a textured quad), rotated text (the run's coverage mask on a rotated quad), text laid out glyph by glyph. |
| `png_out` | Deterministic PNG (no time or text chunks), with the pHYs resolution (and tRNS key) R and Circos write. |
| `text` | Embedded fonts registered with BUI under names of their own (DejaVu Sans for R; CMU Bright Bold, CMU Bright Roman and CMU Sans Serif Demi Condensed for Circos), hinted ink extents from BUI's glyph masks, and R's text advances (the font's hmtx and kern tables at the fractional size, each glyph rounded to whole pixels, as R's Cairo device measures them). |
| `color` | Straight-alpha RGBA; R's colour names and grey levels; R's `hcl()` and `scales::hue_pal()`; `scales::gradient_n_pal()` (interpolation in CIE Lab); Circos's named colours, Brewer lists and `_aN` transparency. |
| `fmath` | sin, cos, pow, log10, cbrt from the `libm` crate, so that every platform computes the same bits. |
| `gg` | ggplot2 4.0 as the hmftools R scripts use it: `theme_bw()` sizes (points, big points, mm line widths at 300 dpi), continuous position scales (limits, expansion, `extended` breaks of the expanded range, log10 breaks, minor breaks, R's `format()` labels), the gtable layout (margins, title, axis titles and labels, panel, guide box), key legends and colour bars, and geoms (rect with R's unanti-aliased fill, point with shape 19, segment and line with R's dash pattern, area, label, stat_bin's bins). `gg::rfmt` holds R's number formatting and the break algorithms. |
| `pdf` | A PDF writer (`Document`: pages, content `Stream`s with PDFBox's operators and number format, TrueType fonts embedded as subsets in a Type0/Identity-H font with a ToUnicode map and the metrics PDFBox computes for a `PDType0Font`, PNG images as image XObjects with soft masks as PDFBox's LosslessFactory makes them (8-bit RGB or grey colour data passed through without recompression), outline, document info), and `pdf::table`, a port of boxable 1.8.2's table layout (tokenizer, line wrapping, broken long words, cell padding, borders, fills, header rows repeated on new pages, rows split across pages, the text-mode and offset handling of its optimized content stream) in the same float arithmetic, so that a Java report's tables come out at the same places. |
| `circos` | A Circos 0.69 engine driven by Circos configuration files: karyotypes (hg19 and hg38 built in, or a file), ideograms (spacing, radius, thickness, fill, stroke, parallel labels), `<highlights>`, plot tracks of type scatter (circle, square, rectangle and triangle glyphs), histogram (joined bins, `extend_bin`), heatmap, highlight, text (Circos's label snuggling) and connector, `<axes>` (`spacing` and `position` lists), `<backgrounds>`, `<rules>` (conditions, `eval()`, `conf()`; on tracks, highlights and links), radii with unit arithmetic (`0.5r-2p`), ideogram `label_format` string expressions, and `<links>` (bezier curves and ribbons with bezier radius, purity and crest). Geometry follows Circos.pm line by line. |

| `gviz` | Gviz's IdeogramTrack as plotTracks draws it (R's Gviz 1.54): bands in Gviz's colours split at the bevel points, rounded caps, centromere triangles, outlines, the range box and the name, at R's png() 72 dpi sizes. |
| `compose` | A `Canvas` (RGB) on which rendered figures are placed unscaled or scaled with bilinear filtering, as cowplot's draw_image and plot_grid compose images (LINX stacks the Circos image, fusion panels and chromosome strip). |
| `gdtext` | GD's string box (libgd 2.3's stringFT bounds) computed from a TrueType font's tables, from which Circos sizes and places labels. |
| `rpdf` | R's pdf() device text metrics for "sans" (Helvetica AFM widths with kerning at whole point sizes): cowplot sizes the ggplot2 layouts it composes on a pdf(NULL) device. |
| `rdev` | R's graphics devices as grid measures text on them, for both faces of "sans": `Device::Png { dpi }` (Cairo with DejaVu Sans Book: hinted advances, fi/fl ligatures, FreeType's synthetic bold), `Device::Pdf` (points; Helvetica and Helvetica-Bold AFM, sizes rounded as R rounds them, kerning as R applies it: pairs naming "space" and each first character's last pair are skipped) and `Device::PdfAt { dpi }` (PDF metrics in pixels, for layouts R sizes on `pdf(NULL)` and draws on png()). String widths, M ascent, per-character metrics, font descent, line spacing, unit conversions. Tables in `rdev_afm.rs` (generated by `scripts/gen_afm_tables.py` from R's AFM files and ISOLatin1 encoding). |
| `gtable` | The gtable package: a table of absolute and null widths and heights with placed grobs (t, l, b, r, z by gtable's rule, clip, name, an optional viewport), adding rows, columns, padding and spacing, subsetting, null units sharing the rest, drawing in z order with clipping, and a dump of every element's device box for comparison with R's viewports. |
| `ggp` | ggplot2's table model for scripts that are matched element for element (ggplot2 3.5.1 or 4.0.3, `ggp::Version`): theme elements, the text model of grid (GEText justification per line, titleGrob margins and ggplot2's descent rule), axes (ticks, labels, guide_axis_logticks), legends (key glyphs, colour bars, guide boxes, legends measured on one device and drawn on another), strips, facet_null and facet_grid (ggh4x's facet_grid2 strip placement and padding), ggplot_gtable's 18 x 15 plot table, geoms (points of shapes 15/19/21, tiles, text, paths with R's dash patterns, lines, ribbons, label keys), scales 1.4.0's labels (`number`, `label_percent`, `precision`, R's `format`) and breaks (`breaks_width`, `breaks_log`, `minor_breaks_log`), and plotmath's `atop(a, atop(b))`. |
| `patchwork` | patchwork 1.2 / 1.3 layouts of ggplot2 tables: designs (`"AABB\nAACC"`), columns and rows with relative heights and widths, nested patchworks, `wrap_elements`, `free()` labels and panels (and `merge()`, which `free()` of a patchwork ends with), strips added to the panel area, aligned panel areas and borders, `plot_annotation` title, subtitle and margins. |
| `pdf_figure` | (features `raster` and `pdf`) display lists written as R's pdf() device writes them: Type 1 Helvetica and Helvetica-Bold (not embedded) with R's encoding differences, text at the rounded size with kerning in TJ arrays, coordinates with 2 decimals, rectangles as `re`, circles as Béziers, alpha by ExtGState, clipping, one page per figure, no dates. |

The plotting code of a tool builds display lists through `gg`, `ggp`, `circos` and `gviz` and never calls BUI.

`text` has `Font::SansBold` (DejaVu Sans with BUI's synthetic-bold flag, the advances unchanged as in FreeType's
emboldening) and maps U+FB01 and U+FB02 (the fi and fl ligatures Pango forms); `color` has R's colour names
mediumseagreen, indianred, lightgrey/lightgray, darkgrey/darkgray and lightcyan2. Neither changes what existing
code draws.

Details of the ggplot2 model that matter for matching R (each checked against R's own viewports or PDF):
- `patchwork::merge` is what `free()` ends with for a patchwork (merge.patchwork): the patchwork nested once more,
  in a 1 x 1 patchwork of its own. It draws the same, but without it the table tree has one level less than R's.
- `patchwork::build_collect` is `build` with `plot_layout(axes = , axis_titles = )` (patchwork 1.3.2's
  collect_axes and collect_axis_titles with merge, after the table is built): of identical axes in a run of patches
  only the outermost is kept and the others' rows or columns emptied (delete_grobs, retrofit_rows / _cols); of
  identical axis titles one is kept, a y title stretched over its run's rows. Grobs are compared by `Placed::key`
  (what R's `hash(unname_grob(x))` compares): the plotting code gives its axes and titles keys that say what they
  draw (`Gtable::set_keys`), and simplify's row and column tables take the keys of what they hold. `build` itself is
  unchanged (`build_collect` without collecting gives the same table). Added for qsee's report (upstream 883ab11bba).
- `patchwork::wrap_full` is `wrap_elements(full = grob)`: the grob over the whole patch, margins included.
- `ggp::geom::key_label` and geom_label boxes are centred on their point (ggplot2 4.0.3's labelGrob raises only the
  text, by `(1 - vjust) * descent`); the key's border width is the layer's linewidth (`KeyGlyph::Label`).
- `ggp::guides::LegendTheme::key_background` (legend.key) is drawn in each key's cell before its glyph, as R does:
  invisible on white, but it covers the edge of the key above.
- `ggp::dash_pattern_on` gives R's dashes on the context's device: Cairo's png() scales the digits by max(lwd, 1),
  pdf() by lwd (`dash_pattern` is Cairo's).
- `ggp::record_drawn_boxes` / `take_drawn_boxes` record, for layout tests, the boxes of what is drawn without a
  table cell of its own (geom_label's rounded rectangles: R's "rrvp" viewports).
- `pdf_figure` writes an axis-aligned rectangle outline as `x y w h re` with each number rounded to 2 decimals, as
  R's pdf() writes rectGrobs, so that the far edges round as R's do.
- `raster`: a clip of zero width or height (the panel of an empty ggplot2 sub-table) shows nothing, as with Cairo
  (it was widened to one pixel). Only `gtable`'s clips can be empty; PURPLE's nine plots are unchanged.

Known differences from R's Cairo PNGs (BUI): glyph origins are whole pixels (Cairo positions glyphs to a fraction
of a pixel), so text differs at glyph edges and a rotated line of text drifts by up to 2 px along its length;
clip rectangles are whole pixels (Cairo clips with partial coverage), so a border drawn inside a clipped panel
differs in its outermost pixel line; line and polygon anti-aliasing is BUI's coverage, not Cairo's.

Opt-in Circos behaviours (off by default, on for LINX): `Circos::with_gd_fill_extent` fills highlights as GD
covers polygons (half a pixel more on each side, so 1-pixel bands are as thick as Circos's), and
`Circos::with_gd_text_bounds` sizes labels by GD's string box rather than by their ink.

## Use

```rust
// a Circos configuration file and the data files beside it
let fig = hmf_plot::circos::draw_conf_file(path, hmf_plot::circos::Options::default())?;
fig.write_png("out.png")?;

// a ggplot2-style plot
let mut fig = hmf_plot::gg::figure(4.8, 4.0); // inches at 300 dpi
let plot = hmf_plot::gg::Plot { /* title, x and y views, axis titles, legends, theme switches */ };
plot.draw(&mut fig, 0.0, 0.0, 1440.0, 1200.0, None, &mut |f, panel| {
    hmf_plot::gg::rect(f, panel, xmin, xmax, ymin, ymax, Some(fill), Some((outline, 0.5)));
});

// a PDF page (feature "pdf"): a content stream with PDFBox's operators, and a boxable table
use hmf_plot::pdf::{Document, PageSize, Rgb, table::Table};
let mut doc = Document::new();
let font = doc.add_font(std::fs::read("font.ttf")?)?;
let page = doc.add_page(PageSize::a4());
let mut s = doc.open_stream(page);
s.begin_text();
s.set_font(font, 10.0);
s.set_non_stroking_color(Rgb::BLACK);
s.new_line_at_offset(30.0, 800.0);
s.show_text(&mut doc, font, "Title")?;
s.end_text();
doc.close_stream(s);
let mut table = Table::new(&doc, 780.0, 741.89, 40.0, 535.0, 30.0, page); // y, y on new pages, bottom, width, x
table.create_row(15.0).create_cell(100.0, "a cell", font)?;
let y_below = table.draw(&mut doc)?; // new pages as needed
doc.save("out.pdf")?;
```

Examples: `cargo run --release --example circos -- <conf> <out.png>` draws a Circos configuration;
`--example side_by_side -- <left.png> <right.png> <out.png> [x y w h zoom]` puts two images side by side.

## Determinism

A PDF has no dates or random ids: the same drawing gives the same bytes (ORANGE's HG008 report is
byte-identical on macOS arm64 and Linux arm64). The same display list gives the same PNG bytes on every run and
on every platform: PURPLE's nine HG008 plots
are byte-identical on macOS arm64, Linux arm64 and Linux amd64. Three things make it so: the fonts are embedded (BUI never falls back to system fonts for
the characters used), the transcendental functions come from `libm`, and the C code is compiled without fused
multiply-adds (`.cargo/config.toml` sets `CFLAGS=-ffp-contract=off`; Apple clang and gcc otherwise fuse
differently inside BUI's polygon coverage, which moved 8 pixels of a Circos image by one unit).

## Building

BUI is a git dependency pinned by `rev` (commit 1af5d5084d92e8ed30c637076bc964c3769386ca of 2026-11-27):
`blenderui-sys` with default features off but `freetype`, which builds the C kit (kit, gpu, core) with CMake.
So the build needs `cmake` and a C/C++ compiler, and network access to github.com (or a vendored copy). The
kit builds in 15-25 s with 4 jobs. The binary then links libc, libm and libgcc_s (and libstdc++ when
referenced).

Cargo prints `error: invalid character '{' in package name: '{{name}}'` once when it resolves the git
dependency: it comes from a template inside BUI's repository and is harmless.

```
CARGO_BUILD_JOBS=4 cargo build --release
CARGO_BUILD_JOBS=4 cargo test --release                 # the raster side (28 tests)
CARGO_BUILD_JOBS=4 cargo test --release --all-features  # and the PDF module (8 more)
cargo test --release --no-default-features --features pdf   # the PDF side alone, without BUI
```

The raster tests (38) check the transformations against what R computes (ggplotGrob widths and heights in points,
breaks, labels, colour ramps, text widths and heights of R's Cairo device, stat_bin bins), the Circos engine's
display lists on small synthetic configurations (angles, radii, colours, rules, joined histograms, label
snuggling, ribbons), and LINX's additions (unit arithmetic, label expressions, glyphs, connectors, axis divisions,
Gviz colours and caps, composition, R's pdf metrics against values from R). The PDF tests check the font metrics
against PDFBox's rules (DejaVu Sans), boxable's line wrapping and page breaks, keeping rows together, the PNG
embedding (pass-through, alpha, tRNS) and the file structure, and that the output is deterministic.

## Licences

BUI: GPL-3.0 (kit sources GPL-2.0-or-later). pdf-writer, subsetter: MIT OR Apache-2.0; ttf-parser: MIT OR
Apache-2.0; miniz_oxide: MIT OR Zlib OR Apache-2.0. DejaVu Sans: Bitstream Vera licence with public-domain changes
(`fonts/LICENSE.DejaVu`). CMU fonts (Bright Bold, Bright Roman, Sans Serif Demi Condensed): X11 licence
(`fonts/LICENSE.CMU`). `src/rpdf_helvetica.rs` and `src/rdev_afm.rs`: metrics from Adobe's Helvetica (and Helvetica-Bold) AFM as R distributes them (Adobe's
notice kept in the file). Circos's karyotype files
(`resources/circos/`): GPL, from Circos 0.69. hmftools is GPL-3.0.

## Upstream

Draws what the R, Circos and matplotlib code of hmftools commit `ca04e6c39d` (2026-11-20) draws; BUI is pinned at `1af5d5084d`. Changes upstream has made since, and what was decided about each, are in
[`../UPSTREAM.md`](../UPSTREAM.md).
