# hmf-plot-rs

Native drawing for the hmftools Rust ports (crate `hmf_plot`): the plots that the Java tools draw with R
(ggplot2) and Circos, drawn in Rust without R, Perl or Circos, and PDF documents. PURPLE's nine plots are the
first user of the raster side (`purple-rs/src/plots.rs`), ORANGE's report the first user of the PDF side
(`orange-rs/src/report/`); LINX's visualiser, CUPPA's and qsee's plots are meant to follow.

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
| `text` | Embedded fonts registered with BUI under names of their own (DejaVu Sans for R; CMU Bright Bold and CMU Sans Serif Demi Condensed for Circos), hinted ink extents from BUI's glyph masks, and R's text advances (the font's hmtx and kern tables at the fractional size, each glyph rounded to whole pixels, as R's Cairo device measures them). |
| `color` | Straight-alpha RGBA; R's colour names and grey levels; R's `hcl()` and `scales::hue_pal()`; `scales::gradient_n_pal()` (interpolation in CIE Lab); Circos's named colours, Brewer lists and `_aN` transparency. |
| `fmath` | sin, cos, pow, log10, cbrt from the `libm` crate, so that every platform computes the same bits. |
| `gg` | ggplot2 4.0 as the hmftools R scripts use it: `theme_bw()` sizes (points, big points, mm line widths at 300 dpi), continuous position scales (limits, expansion, `extended` breaks of the expanded range, log10 breaks, minor breaks, R's `format()` labels), the gtable layout (margins, title, axis titles and labels, panel, guide box), key legends and colour bars, and geoms (rect with R's unanti-aliased fill, point with shape 19, segment and line with R's dash pattern, area, label, stat_bin's bins). `gg::rfmt` holds R's number formatting and the break algorithms. |
| `pdf` | A PDF writer (`Document`: pages, content `Stream`s with PDFBox's operators and number format, TrueType fonts embedded as subsets in a Type0/Identity-H font with a ToUnicode map and the metrics PDFBox computes for a `PDType0Font`, PNG images as image XObjects with soft masks as PDFBox's LosslessFactory makes them (8-bit RGB or grey colour data passed through without recompression), outline, document info), and `pdf::table`, a port of boxable 1.8.2's table layout (tokenizer, line wrapping, broken long words, cell padding, borders, fills, header rows repeated on new pages, rows split across pages, the text-mode and offset handling of its optimized content stream) in the same float arithmetic, so that a Java report's tables come out at the same places. |
| `circos` | A Circos 0.69 engine driven by Circos configuration files: karyotypes (hg19 and hg38 built in, or a file), ideograms (spacing, radius, thickness, fill, stroke, parallel labels), `<highlights>`, plot tracks of type scatter, histogram (joined bins), heatmap, highlight and text (Circos's label snuggling), `<axes>`, `<backgrounds>`, `<rules>` (conditions and `eval()`), and `<links>` (bezier curves and ribbons with bezier radius, purity and crest). Geometry follows Circos.pm line by line. |

The plotting code of a tool builds display lists through `gg` and `circos` and never calls BUI.

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

The raster tests (28) check the transformations against what R computes (ggplotGrob widths and heights in points,
breaks, labels, colour ramps, text widths and heights of R's Cairo device, stat_bin bins) and the Circos
engine's display lists on small synthetic configurations (angles, radii, colours, rules, joined histograms,
label snuggling, ribbons). The PDF tests check the font metrics against PDFBox's rules (DejaVu Sans), boxable's
line wrapping and page breaks, keeping rows together, the PNG embedding (pass-through, alpha, tRNS) and the file
structure, and that the output is deterministic.

## Licences

BUI: GPL-3.0 (kit sources GPL-2.0-or-later). pdf-writer, subsetter: MIT OR Apache-2.0; ttf-parser: MIT OR
Apache-2.0; miniz_oxide: MIT OR Zlib OR Apache-2.0. DejaVu Sans: Bitstream Vera licence with public-domain changes
(`fonts/LICENSE.DejaVu`). CMU fonts: X11 licence (`fonts/LICENSE.CMU`). Circos's karyotype files
(`resources/circos/`): GPL, from Circos 0.69. hmftools is GPL-3.0.
