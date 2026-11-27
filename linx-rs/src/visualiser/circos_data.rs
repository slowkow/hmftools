//! CircosData (the plot's data scaled to Circos positions), CircosConfigWriter (track radii and the
//! configuration from `cluster.template`), CircosDataWriter (the data files) and FusionDataWriter.

use anyhow::Result;
use hmf_common::fmt::{decimal_format, java_double_to_string as jd, java_fixed};
use hmf_common::java_hash::JavaHashSet;

use super::data::*;
use super::files::*;
use super::genome::*;
use super::VisOptions;

pub const CLUSTER_TEMPLATE: &str = include_str!("../../resources/visualisation/cluster.template");

pub const COPY_NUMBER_BASELINE: i32 = 0;
pub const COPY_NUMBER_LOSS_MIN: i32 = -2;
pub const MINOR_ALLELE_COPY_NUMBER_LOSS_MIN: i32 = -1;
pub const AMBER_BAF_MIN: i32 = 0;
pub const AMBER_BAF_MAX: i32 = 1;
pub const COBALT_RATIO_MIN: i32 = 0;
pub const COBALT_RATIO_MAX: i32 = 2;

/// A plot's data after scaling (CircosData).
pub struct CircosData {
    pub exons: Vec<VisGeneExon>,
    pub sv_data: Vec<VisSvData>,
    pub genes: Vec<Gene>,
    pub segments: Vec<VisSegment>,
    pub centromere_sites: Vec<Region>,
    pub line_elements: Vec<Region>,
    pub fragile_sites: Vec<Region>,
    pub copy_numbers: Vec<VisCopyNumber>,
    pub disrupted_gene_regions: Vec<Region>,
    pub connectors: Vec<Connector>,
    pub amber_bafs: Vec<(String, i32, f64)>,
    pub cobalt_ratios: Vec<(String, i32, f64)>,
    pub purple_segments: Vec<PurpleSegment>,
    pub unadjusted_links: Vec<VisSvData>,
    pub unadjusted_copy_numbers: Vec<VisCopyNumber>,
    pub chromosome_ranges: Vec<Region>,
    pub contig_lengths: Vec<GenomePosition>,
    pub upstream_genes: Vec<String>,
    pub downstream_genes: Vec<String>,
    pub copy_number_max: f64,
    pub copy_number_tracks_max: i32,
    pub minor_allele_copy_number_max: f64,
    pub minor_allele_copy_number_tracks_max: i32,
    pub sv_tracks_max: i32,
    pub label_size: f64,
    pub gene_label_size: f64,
    pub max_frame: i32,
    pub gene_relative_size: f64,
}

pub struct CircosInputs<'a> {
    pub segments: Vec<VisSegment>,
    pub links: Vec<VisSvData>,
    pub copy_numbers: Vec<VisCopyNumber>,
    pub exons: &'a [VisGeneExon],
    pub fusions: &'a [VisFusion],
    pub amber: Vec<(String, i32, f64)>,
    pub cobalt: Vec<(String, i32, f64)>,
    pub purple: Vec<PurpleSegment>,
    pub show_simple_sv_segments: bool,
    pub show_fragile_sites: bool,
    pub show_line_elements: bool,
}

fn ceil_round(x: f64) -> i32 {
    java_round(x.ceil()) as i32
}

impl CircosData {
    pub fn new(opts: &VisOptions, sites: &KnownSites, inp: CircosInputs) -> Result<CircosData> {
        let mut upstream_genes: Vec<String> = Vec::new();
        let mut downstream_genes: Vec<String> = Vec::new();
        for f in inp.fusions {
            if !upstream_genes.contains(&f.gene_name_up) {
                upstream_genes.push(f.gene_name_up.clone());
            }
            if !downstream_genes.contains(&f.gene_name_down) {
                downstream_genes.push(f.gene_name_down.clone());
            }
        }
        let unadjusted_disrupted = disrupted_gene_regions(inp.fusions, inp.exons);
        let unadjusted_genes = unique_genes(genes(inp.exons));
        let unadjusted_gene_exons = gene_exons(&unadjusted_genes, inp.exons);

        let mut positions: Vec<GenomePosition> = Vec::new();
        positions.extend(link_positions(&inp.links));
        positions.extend(all_positions(&inp.segments));
        if opts.interpolate_copy_number_positions {
            positions.extend(min_max_positions(&inp.copy_numbers));
        } else {
            positions.extend(all_positions(&inp.copy_numbers));
        }
        if !opts.interpolate_exon_positions {
            positions.extend(all_positions(&unadjusted_gene_exons));
        }
        let mut positions = distinct_positions(positions);
        sort_positions(&mut positions);

        let span = span_positions(&positions);
        let unadjusted_centromeres = limit_highlights(&sites.centromeres, &span);
        let unadjusted_fragile = if inp.show_fragile_sites { limit_highlights(&sites.fragile_sites, &span) } else { Vec::new() };
        let unadjusted_line = if inp.show_line_elements { limit_highlights(&sites.line_elements, &span) } else { Vec::new() };

        let scaler = ScalePosition::new(&positions);

        let mut segments = Vec::new();
        for s in &inp.segments {
            let c = scaler.contig(&s.chromosome)?;
            let mut n = s.clone();
            n.start = c.scale(s.start)?;
            n.end = c.scale(s.end)?;
            segments.push(n);
        }
        let mut sv_data = Vec::new();
        for l in &inp.links {
            let mut n = l.clone();
            if l.is_valid_start() {
                n.pos_start = scaler.contig(&l.chr_start)?.scale(l.pos_start)?;
            }
            if l.is_valid_end() {
                n.pos_end = scaler.contig(&l.chr_end)?.scale(l.pos_end)?;
            }
            sv_data.push(n);
        }
        let mut copy_numbers = Vec::new();
        for cn in &inp.copy_numbers {
            let c = scaler.contig(&cn.chromosome)?;
            let mut n = cn.clone();
            n.start = c.scale(cn.start)?;
            n.end = c.scale(cn.end)?;
            copy_numbers.push(n);
        }
        let interp_regions = |rs: &[Region]| -> Result<Vec<Region>> {
            rs.iter()
                .map(|r| {
                    let c = scaler.contig(&r.chromosome)?;
                    Ok(Region::new(&r.chromosome, c.interpolate(r.start), c.interpolate(r.end)))
                })
                .collect()
        };
        let centromere_sites = interp_regions(&unadjusted_centromeres)?;
        let line_elements = interp_regions(&unadjusted_line)?;
        let fragile_sites = interp_regions(&unadjusted_fragile)?;

        // ScalePosition.interpolateGene
        let gene_name_distance = (1.5 / 360.0) * scaler.total_length as f64;
        let mut genes_scaled = Vec::new();
        for g in &unadjusted_genes {
            let c = scaler.contig(&g.chromosome)?;
            let s = c.interpolate(g.start);
            let e = c.interpolate(g.end);
            let mut name_pos = java_round(if g.strand > 0 { s as f64 - gene_name_distance } else { e as f64 + gene_name_distance }) as i32;
            if name_pos < 1 {
                name_pos = java_round(e as f64 + gene_name_distance) as i32;
            }
            if name_pos > c.length() {
                name_pos = java_round(s as f64 - gene_name_distance) as i32;
            }
            let mut n = g.clone();
            n.start = s;
            n.end = e;
            n.name_position = name_pos as i64;
            genes_scaled.push(n);
        }
        let disrupted_gene_regions = interp_regions(&unadjusted_disrupted)?;
        let mut exons = Vec::new();
        for e in &unadjusted_gene_exons {
            let c = scaler.contig(&e.chromosome)?;
            let mut n = e.clone();
            n.exon_start = c.interpolate(e.exon_start);
            n.exon_end = c.interpolate(e.exon_end);
            exons.push(n);
        }

        let interp_points = |pts: Vec<(String, i32, f64)>| -> Result<Vec<(String, i32, f64)>> {
            pts.into_iter().map(|(chr, p, v)| Ok((chr.clone(), scaler.contig(&chr)?.interpolate(p), v))).collect()
        };
        let amber_bafs = interp_points(downsample_points(&inp.amber, &positions))?;
        let cobalt_ratios = interp_points(center_cobalt_positions(&downsample_points(&inp.cobalt, &positions), &positions))?;
        let mut purple_segments = Vec::new();
        for p in &inp.purple {
            let c = scaler.contig(&p.chromosome)?;
            let mut n = p.clone();
            n.pos_start = c.interpolate(p.max_start);
            n.pos_end = c.interpolate(p.pos_end);
            purple_segments.push(n);
        }

        let sv_tracks_max = segments.iter().map(|s| s.track).max().unwrap_or(0) + 1;
        let copy_number_max = copy_numbers.iter().map(|c| c.copy_number).fold(None, |m: Option<f64>, v| Some(m.map_or(v, |m| m.max(v)))).unwrap_or(0.0);
        let copy_number_tracks_max = i32::max(2, ceil_round(copy_number_max - 2.0));
        let minor_allele_copy_number_max =
            copy_numbers.iter().map(|c| c.minor_allele_copy_number()).fold(None, |m: Option<f64>, v| Some(m.map_or(v, |m| m.max(v)))).unwrap_or(0.0);
        let minor_allele_copy_number_tracks_max = i32::max(1, ceil_round(minor_allele_copy_number_max - 1.0));

        let connectors = create_connectors(inp.show_simple_sv_segments, &segments, &sv_data);
        let untruncated = copy_numbers.iter().filter(|c| !c.truncated).count() as i64;
        let label_size = opts.label_size(untruncated) as f64;
        let max_gene_chars = genes_scaled.iter().map(|g| g.name.encode_utf16().count() as i32).max().unwrap_or(0);
        let gene_label_size = if max_gene_chars > opts.max_gene_characters {
            0.9 * opts.max_gene_characters as f64 / max_gene_chars as f64 * label_size
        } else {
            label_size
        };
        let max_frame = segments.iter().map(|s| s.frame).max().unwrap_or(0);

        Ok(CircosData {
            exons,
            sv_data,
            genes: genes_scaled,
            segments,
            centromere_sites,
            line_elements,
            fragile_sites,
            copy_numbers,
            disrupted_gene_regions,
            connectors,
            amber_bafs,
            cobalt_ratios,
            purple_segments,
            unadjusted_links: inp.links,
            unadjusted_copy_numbers: inp.copy_numbers,
            chromosome_ranges: span,
            contig_lengths: scaler.contig_lengths.clone(),
            upstream_genes,
            downstream_genes,
            copy_number_max,
            copy_number_tracks_max,
            minor_allele_copy_number_max,
            minor_allele_copy_number_tracks_max,
            sv_tracks_max,
            label_size,
            gene_label_size,
            max_frame,
            gene_relative_size: opts.gene_relative_size,
        })
    }

    pub fn display_genes(&self) -> bool {
        !self.exons.is_empty() && self.gene_relative_size - 0.0 > 1e-10
    }

    pub fn total_contig_length(&self) -> i32 {
        self.contig_lengths.iter().map(|c| c.position).sum()
    }
}

// ---------------------------------------------------------------------------------------------------------
// CircosConfigWriter

pub struct ConfigWriter {
    exon_outer: f64,
    exon_inner: f64,
    gene_outer: f64,
    gene_inner: f64,
    cobalt_outer: f64,
    cobalt_inner: f64,
    amber_outer: f64,
    amber_inner: f64,
    segment_outer: f64,
    segment_inner: f64,
    cn_outer: f64,
    cn_middle: f64,
    cn_inner: f64,
    map_outer: f64,
    map_middle: f64,
    map_inner: f64,
    innermost: f64,
    label_size: f64,
    label_position: f64,
}

impl ConfigWriter {
    pub fn new(data: &CircosData, opts: &VisOptions, amber: bool, cobalt: bool) -> ConfigWriter {
        let gap = opts.gap_radius;
        let display_genes = data.display_genes();
        let gene_rel = if display_genes { opts.gene_relative_size } else { 0.0 };
        let seg_rel = opts.segment_relative_size;
        let cn_rel = opts.copy_number_relative_size;
        let total_rel = gene_rel + seg_rel + cn_rel;
        let mut gaps = 3;
        if display_genes {
            gaps += 2;
        }
        if amber {
            gaps += 1;
        }
        if cobalt {
            gaps += 1;
        }
        let total_space = 1.0 - gaps as f64 * gap - opts.inner_radius;
        let purple_space = cn_rel / total_rel * total_space;
        let cna_gain_rel = i32::max(2, ceil_round(data.copy_number_max - 2.0));
        let cna_loss_rel = 2;
        let map_gain_rel = i32::max(1, ceil_round(data.minor_allele_copy_number_max - 1.0));
        let map_loss_rel = 1;
        let amber_rel = if amber { 2 } else { 0 };
        let cobalt_rel = if cobalt { 2 } else { 0 };
        let purple_rel = cna_gain_rel + cna_loss_rel + map_gain_rel + map_loss_rel + amber_rel + cobalt_rel;
        let purple_track = purple_space / purple_rel as f64;

        let mut w = ConfigWriter {
            exon_outer: 0.0,
            exon_inner: 0.0,
            gene_outer: 0.0,
            gene_inner: 0.0,
            cobalt_outer: 0.0,
            cobalt_inner: 0.0,
            amber_outer: 0.0,
            amber_inner: 0.0,
            segment_outer: 0.0,
            segment_inner: 0.0,
            cn_outer: 0.0,
            cn_middle: 0.0,
            cn_inner: 0.0,
            map_outer: 0.0,
            map_middle: 0.0,
            map_inner: 0.0,
            innermost: 0.0,
            label_size: data.label_size,
            label_position: 0.0,
        };
        let mut current = 1.0;
        if display_genes {
            w.exon_outer = current - gap - opts.exon_rank_radius;
            w.exon_inner = w.exon_outer - gene_rel / total_rel * total_space;
            current = w.exon_inner;
            let exon_track = w.exon_outer - w.exon_inner;
            let trim = 9.0 / 20.0 * exon_track;
            w.gene_outer = w.exon_outer - trim;
            w.gene_inner = w.exon_inner + trim;
        }
        if cobalt {
            w.cobalt_outer = current - gap;
            w.cobalt_inner = w.cobalt_outer - cobalt_rel as f64 * purple_track;
            current = w.cobalt_inner;
        }
        if amber {
            w.amber_outer = current - gap;
            w.amber_inner = w.amber_outer - amber_rel as f64 * purple_track;
            current = w.amber_inner;
        }
        w.segment_outer = current - gap;
        w.segment_inner = w.segment_outer - seg_rel / total_rel * total_space;
        current = w.segment_inner;
        w.cn_outer = current - gap;
        w.cn_middle = w.cn_outer - cna_gain_rel as f64 * purple_track;
        w.cn_inner = w.cn_middle - cna_loss_rel as f64 * purple_track;
        current = w.cn_inner;
        w.map_outer = current - gap;
        w.map_middle = w.map_outer - map_gain_rel as f64 * purple_track;
        w.map_inner = w.map_middle - map_loss_rel as f64 * purple_track;
        current = w.map_inner;
        w.innermost = current;
        w.label_position = w.cn_middle + purple_track * 0.05;
        w
    }

    pub fn sv_track_relative(&self, track: i32, sv_tracks_max: i32) -> f64 {
        let difference = self.segment_outer - self.segment_inner;
        let single = difference / sv_tracks_max as f64;
        self.segment_inner + track as f64 * single
    }

    /// writeConfig(frame): the template with its SUBSTITUTE_ fields replaced, in Java's order
    pub fn config(&self, data: &CircosData, opts: &VisOptions, file_id: &str, frame: i32) -> String {
        let chromosome_count = data.contig_lengths.len();
        let total = data.total_contig_length();
        let spacing = if chromosome_count > 1 { "0.005r".to_string() } else { format!("{}u", jd(0.005 * total as f64)) };
        let reps: Vec<(&str, String)> = vec![
            ("SUBSTITUTE_CURRENT_FRAME", frame.to_string()),
            ("SUBSTITUTE_IDEOGRAM_RADIUS", jd(opts.outer_radius)),
            ("SUBSTITUTE_IDEOGRAM_SPACING", spacing),
            ("SUBSTITUTE_IDEOGRAM_THICKNESS", "50".to_string()),
            ("SUBSTITUTE_EXON_INNER_RADIUS", jd(self.exon_inner)),
            ("SUBSTITUTE_EXON_OUTER_RADIUS", jd(self.exon_outer)),
            ("SUBSTITUTE_GENE_INNER_RADIUS", jd(self.gene_inner)),
            ("SUBSTITUTE_GENE_OUTER_RADIUS", jd(self.gene_outer)),
            ("SUBSTITUTE_COBALT_INNER_RADIUS", jd(self.cobalt_inner)),
            ("SUBSTITUTE_COBALT_OUTER_RADIUS", jd(self.cobalt_outer)),
            ("SUBSTITUTE_COBALT_RATIO_MAX", COBALT_RATIO_MAX.to_string()),
            ("SUBSTITUTE_COBALT_RATIO_MIN", COBALT_RATIO_MIN.to_string()),
            ("SUBSTITUTE_AMBER_INNER_RADIUS", jd(self.amber_inner)),
            ("SUBSTITUTE_AMBER_OUTER_RADIUS", jd(self.amber_outer)),
            ("SUBSTITUTE_AMBER_BAF_MAX", AMBER_BAF_MAX.to_string()),
            ("SUBSTITUTE_AMBER_BAF_MIN", AMBER_BAF_MIN.to_string()),
            ("SUBSTITUTE_SV_INNER_RADIUS", jd(self.segment_inner)),
            ("SUBSTITUTE_SV_OUTER_RADIUS", jd(self.segment_outer)),
            ("SUBSTITUTE_MINOR_ALLELE_COPY_NUMBER_INNER_RADIUS", jd(self.map_inner)),
            ("SUBSTITUTE_MINOR_ALLELE_COPY_NUMBER_OUTER_RADIUS", jd(self.map_outer)),
            ("SUBSTITUTE_MINOR_ALLELE_COPY_NUMBER_MIDDLE_RADIUS", jd(self.map_middle)),
            ("SUBSTITUTE_MINOR_ALLELE_COPY_NUMBER_LOSS_MIN", MINOR_ALLELE_COPY_NUMBER_LOSS_MIN.to_string()),
            ("SUBSTITUTE_MINOR_ALLELE_COPY_NUMBER_GAIN_MAX", data.minor_allele_copy_number_tracks_max.to_string()),
            ("SUBSTITUTE_MINOR_ALLELE_COPY_NUMBER_GAIN_SPACING", jd(1.0 / data.minor_allele_copy_number_tracks_max as f64)),
            ("SUBSTITUTE_COPY_NUMBER_INNER_RADIUS", jd(self.cn_inner)),
            ("SUBSTITUTE_COPY_NUMBER_OUTER_RADIUS", jd(self.cn_outer)),
            ("SUBSTITUTE_COPY_NUMBER_MIDDLE_RADIUS", jd(self.cn_middle)),
            ("SUBSTITUTE_COPY_NUMBER_BASELINE", COPY_NUMBER_BASELINE.to_string()),
            ("SUBSTITUTE_COPY_NUMBER_LOSS_MIN", COPY_NUMBER_LOSS_MIN.to_string()),
            ("SUBSTITUTE_COPY_NUMBER_GAIN_MAX", data.copy_number_tracks_max.to_string()),
            ("SUBSTITUTE_COPY_NUMBER_GAIN_AXIS_POSITION", cna_axis_positions(data.copy_number_tracks_max)),
            ("SUBSTITUTE_COPY_NUMBER_DISTANCE_RADIUS", format!("{}r", jd(self.label_position))),
            ("SUBSTITUTE_INNERMOST_RADIUS", jd(self.innermost)),
            ("SUBSTITUTE_SV_SPACING", jd(1.0 / data.sv_tracks_max as f64)),
            ("SUBSTITUTE_SV_MAX", data.sv_tracks_max.to_string()),
            ("SUBSTITUTE_LABEL_SIZE", jd(self.label_size)),
            ("SUBSTITUTE_SAMPLE", file_id.to_string()),
        ];
        let mut t = CLUSTER_TEMPLATE.to_string();
        for (k, v) in reps {
            t = t.replace(k, &v);
        }
        t
    }
}

/// CircosConfigWriter.cnaAxisPositions
pub fn cna_axis_positions(max_tracks: i32) -> String {
    let rel = 1.0 / max_tracks as f64;
    let s = |i: i32| format!("{}r", jd(java_round(i as f64 * rel * 10000.0) as f64 / 10000.0));
    let mut parts = Vec::new();
    for i in 1..=i32::min(7, max_tracks) {
        parts.push(s(i));
    }
    let mut i = 8;
    while i <= max_tracks {
        parts.push(s(i));
        i += 10;
    }
    parts.join(",")
}

// ---------------------------------------------------------------------------------------------------------
// CircosDataWriter

const SINGLE_BLUE: &str = "(107,174,214)";
const SINGLE_RED: &str = "(214,144,107)";
const SINGLE_GREEN: &str = "(107,214,148)";
const SINGLE_YELLOW: &str = "(214,210,107)";
const MAX_CONTIG_LENGTH_TO_DISPLAY_EXON_RANK: i32 = 100000;
const MIN_KARYOTYPE_LENGTH: i32 = 10;

/// CircosDataWriter.shorthand
pub fn shorthand(value: i32) -> String {
    if value < 100 {
        value.to_string()
    } else if value < 99_950 {
        format!("{}k", java_fixed(value as f64 / 1000.0, 1))
    } else {
        format!("{}m", java_fixed(value as f64 / 1_000_000.0, 1))
    }
}

/// DecimalFormat("#,###") (Locale default: the grouping comma)
fn position_format(v: i32) -> String {
    let s = v.unsigned_abs().to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    if v < 0 { format!("-{out}") } else { out }
}

/// The data files of a plot, by suffix (as Java writes them: each line followed by a newline).
pub struct CircosFiles {
    pub files: Vec<(&'static str, Vec<String>)>,
}

impl CircosFiles {
    pub fn get(&self, suffix: &str) -> Option<&Vec<String>> {
        self.files.iter().find(|(s, _)| *s == suffix).map(|(_, l)| l)
    }
}

pub fn write_circos_data(data: &CircosData, colours: &ColorPicker, cw: &ConfigWriter, opts: &VisOptions) -> CircosFiles {
    let thickness = Thickness::new(opts.min_line_size as f64, opts.max_line_size as f64, &data.connectors);
    let mut gene_colours: std::collections::HashMap<String, &str> = std::collections::HashMap::new();
    for g in &data.genes {
        match g.annotation_type {
            AnnotationType::PSEUDOGENE => {
                gene_colours.insert(g.name.clone(), SINGLE_YELLOW);
            }
            AnnotationType::DRIVER | AnnotationType::DISRUPTION => {
                gene_colours.insert(g.name.clone(), SINGLE_GREEN);
            }
            AnnotationType::FUSION => {
                gene_colours.insert(g.name.clone(), SINGLE_BLUE);
            }
            AnnotationType::EXON_LOST => {}
        }
    }
    for g in &data.upstream_genes {
        gene_colours.insert(g.clone(), SINGLE_BLUE);
    }
    for g in &data.downstream_genes {
        gene_colours.insert(g.clone(), SINGLE_RED);
    }
    let gc = |name: &str| gene_colours.get(name).copied().unwrap_or("null");
    let total = data.total_contig_length();
    let display_genes = data.display_genes();
    let tab = |v: &[String]| v.join("\t");

    let mut exon_lines = Vec::new();
    let mut rank_lines = Vec::new();
    let mut gene_lines = Vec::new();
    let mut gene_name_lines = Vec::new();
    if display_genes {
        for e in &data.exons {
            exon_lines.push(tab(&[circos_contig(&e.chromosome), e.exon_start.to_string(), e.exon_end.to_string(), "1".into(), format!("fill_color={}", gc(&e.gene))]));
        }
        for r in &data.disrupted_gene_regions {
            exon_lines.push(tab(&[circos_contig(&r.chromosome), r.start.to_string(), r.end.to_string(), "1".into(), "fill_color=(255,255,255,0.6)".into()]));
        }
        if total <= MAX_CONTIG_LENGTH_TO_DISPLAY_EXON_RANK {
            for e in &data.exons {
                let p = e.exon_start + (e.exon_end - e.exon_start) / 2;
                rank_lines.push(tab(&[circos_contig(&e.chromosome), p.to_string(), p.to_string(), e.exon_rank.to_string()]));
            }
        }
        for g in &data.genes {
            gene_lines.push(tab(&[circos_contig(&g.chromosome), g.start.to_string(), g.end.to_string(), "1".into(), format!("fill_color={}", gc(&g.name))]));
            gene_name_lines.push(tab(&[
                circos_contig(&g.chromosome),
                g.name_position.to_string(),
                g.name_position.to_string(),
                g.name.clone(),
                format!("label_size={}p", jd(data.gene_label_size)),
            ]));
        }
    }

    let position_lines = position_text(&data.unadjusted_links, &data.sv_data, opts);

    let mut segment_lines = Vec::new();
    for s in &data.segments {
        if s.track > 0 {
            let t = thickness.pixels(s.link_ploidy);
            let r0 = cw.sv_track_relative(s.track, data.sv_tracks_max);
            let r = decimal_format(r0, 3);
            segment_lines.push(tab(&[
                circos_contig(&s.chromosome),
                s.start.to_string(),
                s.end.to_string(),
                s.track.to_string(),
                format!(
                    "fill_{},r0={}r-{}p,r1={}r+{}p,frame={}",
                    colours.transparent_colour(s.cluster_id, s.chain_id),
                    r,
                    jd(t / 2.0),
                    r,
                    jd(t / 2.0),
                    s.frame
                ),
            ]));
        }
    }

    let mut sorted_contigs = data.contig_lengths.clone();
    sort_positions(&mut sorted_contigs);
    let karyotype_lines: Vec<String> = sorted_contigs
        .iter()
        .map(|c| {
            let name = human_chromosome_name(&c.chromosome).unwrap_or(&c.chromosome);
            format!("chr - {} {} 1 {} chr{}", circos_contig(&c.chromosome), name, i32::max(MIN_KARYOTYPE_LENGTH, c.position), name)
        })
        .collect();

    let connector_lines: Vec<String> = data
        .connectors
        .iter()
        .map(|c| {
            let r1 = cw.sv_track_relative(c.track, data.sv_tracks_max);
            tab(&[
                circos_contig(&c.chromosome),
                c.position.to_string(),
                c.position.to_string(),
                format!(
                    "r1={}r,{},thickness={},frame={}",
                    decimal_format(r1, 3),
                    colours.transparent_colour(c.cluster_id, c.chain_id),
                    jd(thickness.pixels(c.ploidy)),
                    c.frame
                ),
            ])
        })
        .collect();

    let mut link_lines = Vec::new();
    for l in &data.sv_data {
        if !human_chromosome_contains(&l.chr_start) || !human_chromosome_contains(&l.chr_end) {
            continue;
        }
        link_lines.push(tab(&[
            circos_contig(&l.chr_start),
            l.pos_start.to_string(),
            l.pos_start.to_string(),
            circos_contig(&l.chr_end),
            l.pos_end.to_string(),
            l.pos_end.to_string(),
            format!("{},thickness={},frame={}", colours.transparent_colour(l.cluster_id, l.chain_id), jd(thickness.pixels(l.jcn)), l.frame),
        ]));
    }

    // createScatter
    let glyph_size = opts.glyph_size;
    let glyph_inner = (opts.glyph_size as f64 * 14.0 / 20.0).floor() as i32;
    let mut scatter_lines = Vec::new();
    for s in &data.segments {
        if s.track == 0 {
            continue;
        }
        let colour = colours.transparent_colour(s.cluster_id, s.chain_id);
        for is_start in [true, false] {
            let loc = if is_start { s.start } else { s.end };
            let terminal = if is_start { s.terminal_start } else { s.terminal_end };
            let glyph = if terminal != SegmentTerminal::NONE {
                "square"
            } else {
                let fb = find_start_link(&s.chromosome, loc, &data.sv_data).is_some_and(|l| l.info_start == "FOLDBACK")
                    || find_end_link(&s.chromosome, loc, &data.sv_data).is_some_and(|l| l.info_end == "FOLDBACK");
                if fb { "triangle" } else { "circle" }
            };
            let entry = |c: &str, size: i32, frame: i32| {
                tab(&[
                    circos_contig(&s.chromosome),
                    loc.to_string(),
                    loc.to_string(),
                    s.track.to_string(),
                    format!("{c},glyph={glyph},glyph_size={size},frame={frame}"),
                ])
            };
            scatter_lines.push(entry(&colour, glyph_size, s.frame));
            if terminal == SegmentTerminal::CENTROMERE {
                scatter_lines.push(entry("color=white", glyph_inner, 0));
            }
        }
    }
    let mut sgl_lines = Vec::new();
    for l in &data.sv_data {
        if l.is_valid_start() && !l.is_valid_end() {
            let colour = colours.transparent_colour(l.cluster_id, l.chain_id);
            for (c, size) in [(colour.as_str(), glyph_size), ("color=white", glyph_inner)] {
                sgl_lines.push(tab(&[
                    circos_contig(&l.chr_start),
                    l.pos_start.to_string(),
                    l.pos_start.to_string(),
                    "0".into(),
                    format!("{c},glyph=circle,glyph_size={size}"),
                ]));
            }
        }
    }

    let cna_lines: Vec<String> =
        data.copy_numbers.iter().map(|c| tab(&[circos_contig(&c.chromosome), c.start.to_string(), c.end.to_string(), jd(c.copy_number - 2.0)])).collect();
    let map_lines: Vec<String> = data
        .copy_numbers
        .iter()
        .map(|c| tab(&[circos_contig(&c.chromosome), c.start.to_string(), c.end.to_string(), jd(c.minor_allele_copy_number() - 1.0)]))
        .collect();
    let highlights = |rs: &[Region]| -> Vec<String> { rs.iter().map(|r| tab(&[circos_contig(&r.chromosome), r.start.to_string(), r.end.to_string()])).collect() };

    let mut distance_lines = Vec::new();
    let untruncated = data.copy_numbers.iter().filter(|c| !c.truncated).count() as i32;
    if untruncated <= opts.max_distance_labels {
        for (i, unadjusted) in data.unadjusted_copy_numbers.iter().enumerate() {
            let adjusted = &data.copy_numbers[i];
            if !adjusted.truncated {
                distance_lines.push(tab(&[
                    circos_contig(&adjusted.chromosome),
                    adjusted.start.to_string(),
                    adjusted.end.to_string(),
                    shorthand(unadjusted.bases()),
                ]));
            }
        }
    }

    let chromosome_lines: Vec<String> =
        data.chromosome_ranges.iter().map(|r| tab(&[r.chromosome.clone(), r.start.to_string(), r.end.to_string(), hex_contig_colour(&r.chromosome)])).collect();
    let amber_lines: Vec<String> = data.amber_bafs.iter().map(|(c, p, v)| tab(&[circos_contig(c), p.to_string(), p.to_string(), jd(*v)])).collect();
    let amber_pcf_lines: Vec<String> = data
        .purple_segments
        .iter()
        .map(|s| tab(&[circos_contig(&s.chromosome), s.pos_start.to_string(), s.pos_end.to_string(), jd(f64::min(1.0 - s.observed_baf, s.observed_baf))]))
        .collect();
    let cobalt_lines: Vec<String> = data.cobalt_ratios.iter().map(|(c, p, v)| tab(&[circos_contig(c), p.to_string(), p.to_string(), jd(*v)])).collect();
    let cobalt_pcf_lines: Vec<String> = data
        .purple_segments
        .iter()
        .filter(|s| s.germline_diploid)
        .map(|s| tab(&[circos_contig(&s.chromosome), s.pos_start.to_string(), s.pos_end.to_string(), jd(s.observed_tumor_ratio)]))
        .collect();

    CircosFiles {
        files: vec![
            ("exon.circos", exon_lines),
            ("exon.rank.circos", rank_lines),
            ("gene.circos", gene_lines),
            ("gene.name.circos", gene_name_lines),
            ("position.circos", position_lines),
            ("segment.circos", segment_lines),
            ("karyotype.circos", karyotype_lines),
            ("connector.circos", connector_lines),
            ("link.circos", link_lines),
            ("scatter.circos", scatter_lines),
            ("scatter.sgl.circos", sgl_lines),
            ("cna.circos", cna_lines),
            ("map.circos", map_lines),
            ("centromeres.circos", highlights(&data.centromere_sites)),
            ("fragile.circos", highlights(&data.fragile_sites)),
            ("line_element.circos", highlights(&data.line_elements)),
            ("distance.circos", distance_lines),
            ("chromosome.circos", chromosome_lines),
            ("amber.circos", amber_lines),
            ("amber_pcf.circos", amber_pcf_lines),
            ("cobalt.circos", cobalt_lines),
            ("cobalt_pcf.circos", cobalt_pcf_lines),
        ],
    }
}

/// createPositionText: position labels at least 100k (else 1M, else 10M) apart, sorted as strings
fn position_text(original: &[VisSvData], scaled: &[VisSvData], opts: &VisOptions) -> Vec<String> {
    // AdjustedPositions.create
    let mut positions: Vec<(String, i32, i32, i32)> = Vec::new(); // chromosome, position, unadjusted, svId
    for (o, s) in original.iter().zip(scaled) {
        if s.is_valid_start() {
            positions.push((s.chr_start.clone(), s.pos_start, o.pos_start, o.sv_id));
        }
        if s.is_valid_end() {
            positions.push((s.chr_end.clone(), s.pos_end, o.pos_end, o.sv_id));
        }
    }
    positions.sort_by(|a, b| contig_cmp(&a.0, &b.0).then(a.1.cmp(&b.1)));
    let make = |min_distance: i32, exact: bool| -> Vec<String> {
        let mut result: Vec<String> = Vec::new();
        let mut contigs: Vec<&str> = Vec::new();
        for p in &positions {
            if !contigs.contains(&p.0.as_str()) {
                contigs.push(&p.0);
            }
        }
        for contig in contigs {
            let mut current = 0;
            for p in positions.iter().filter(|p| p.0 == contig) {
                let new_pos = p.2;
                if new_pos - min_distance >= current {
                    let mut label = if exact { position_format(p.2) } else { shorthand(p.2) };
                    if opts.show_sv_id {
                        label += &format!(":{}", p.3);
                    }
                    result.push([circos_contig(contig), p.1.to_string(), p.1.to_string(), label].join("\t"));
                    current = new_pos;
                }
            }
        }
        result.sort();
        result.dedup();
        result
    };
    if opts.exact_position {
        return make(1, true);
    }
    let a = make(100_000, false);
    if (a.len() as i32) < opts.max_position_labels {
        return a;
    }
    let b = make(1_000_000, false);
    if (b.len() as i32) < opts.max_position_labels {
        return b;
    }
    make(10_000_000, false)
}

// ---------------------------------------------------------------------------------------------------------
// FusionDataWriter

#[derive(Clone, Debug)]
pub struct FusedExon {
    pub cluster_id: i32,
    pub fusion: String,
    pub gene: String,
    pub gene_start: i32,
    pub gene_end: i32,
    pub is_upstream: bool,
    pub chromosome: String,
    pub start: i32,
    pub end: i32,
    pub rank: i32,
    pub skipped: bool,
    pub transcript: String,
}

/// FusedExons.convertRegion
pub(crate) fn convert_region(strand: i32, reference: &Region, start: i32, end: i32) -> (i32, i32) {
    let (s, e) = if strand < 0 {
        (reference.end - end, reference.end - start)
    } else {
        (start - reference.start, i32::min(reference.end, end) - reference.start)
    };
    (i32::max(0, s), e)
}

fn overlaps(a: &Region, s: i32, e: i32) -> bool {
    e > a.start && s < a.end
}

/// FusedExons.fusedExons
pub fn fused_exons(f: &VisFusion, exons: &[VisGeneExon]) -> Vec<FusedExon> {
    let fe = fusion_exons(f, exons);
    let up = sorted_upstream_exons(f, &fe);
    let down = sorted_downstream_exons(f, &fe);
    if up.is_empty() || down.is_empty() {
        return Vec::new();
    }
    let first_up = &up[0];
    let up_region = if f.strand_up > 0 {
        Region::new(&first_up.chromosome, first_up.exon_start, f.pos_up)
    } else {
        Region::new(&first_up.chromosome, f.pos_up, first_up.exon_end)
    };
    let conv_up = convert_region(f.strand_up, &up_region, up_region.start, up_region.end);
    let mut result = Vec::new();
    let mut has_up = false;
    for e in &up {
        let c = convert_region(f.strand_up, &up_region, e.exon_start, e.exon_end);
        if overlaps(&up_region, e.exon_start, e.exon_end) && e.chromosome == up_region.chromosome {
            result.push(FusedExon {
                cluster_id: f.cluster_id,
                fusion: f.name(),
                gene: f.gene_name_up.clone(),
                gene_start: conv_up.0,
                gene_end: conv_up.1,
                is_upstream: true,
                chromosome: f.chr_up.clone(),
                start: c.0,
                end: c.1,
                rank: e.exon_rank,
                skipped: e.exon_rank > f.fused_exon_up,
                transcript: f.transcript_up.clone(),
            });
            has_up = true;
        }
    }
    let final_down = &down[down.len() - 1];
    let down_region = if f.strand_down > 0 {
        Region::new(&final_down.chromosome, f.pos_down, final_down.exon_end)
    } else {
        Region::new(&final_down.chromosome, final_down.exon_start, f.pos_down)
    };
    let conv_down = convert_region(f.strand_down, &down_region, down_region.start, down_region.end);
    let intronic_to_exonic = f.region_type_up == "Intronic" && f.region_type_down == "Exonic";
    let mut has_down = false;
    for (i, e) in down.iter().enumerate() {
        let c = convert_region(f.strand_down, &down_region, e.exon_start, e.exon_end);
        if overlaps(&down_region, e.exon_start, e.exon_end) && e.chromosome == down_region.chromosome {
            result.push(FusedExon {
                cluster_id: f.cluster_id,
                fusion: f.name(),
                gene: f.gene_name_down.clone(),
                gene_start: conv_down.0 + conv_up.1,
                gene_end: conv_down.1 + conv_up.1,
                is_upstream: false,
                chromosome: f.chr_down.clone(),
                start: c.0 + conv_up.1,
                end: c.1 + conv_up.1,
                rank: e.exon_rank,
                skipped: e.exon_rank < f.fused_exon_down || (i == 0 && intronic_to_exonic),
                transcript: f.transcript_down.clone(),
            });
            has_down = true;
        }
    }
    if !has_up || !has_down {
        return Vec::new();
    }
    result
}

/// VisProteinDomains.exonicProteinDomains: UTR domains cut to the exons they overlap
pub fn exonic_protein_domains(domains: &[VisProteinDomain], exons: &[VisGeneExon]) -> Vec<VisProteinDomain> {
    let mut out = Vec::new();
    for d in domains {
        if d.name() == PROTEIN_DOMAIN_UTR {
            let mut pe: Vec<&VisGeneExon> = exons.iter().filter(|e| e.transcript == d.transcript).collect();
            pe.sort_by(|a, b| region_cmp(*a, *b));
            for e in pe {
                if d.chromosome == e.chromosome && e.exon_end > d.start && e.exon_start < d.end {
                    let mut n = d.clone();
                    n.start = i32::max(n.start, e.exon_start);
                    n.end = i32::min(n.end, e.exon_end);
                    out.push(n);
                }
            }
        } else {
            out.push(d.clone());
        }
    }
    out
}

/// FusedProteinDomains.fusedProteinDomains
pub fn fused_protein_domains(f: &VisFusion, fused: &[FusedExon], domains: &[VisProteinDomain]) -> Vec<VisProteinDomain> {
    let mut out = Vec::new();
    if fused.is_empty() {
        return out;
    }
    let first_up = &fused[0];
    let up_len = first_up.gene_end - first_up.gene_start;
    let up_start = if f.strand_up < 0 { f.pos_up } else { f.pos_up - up_len };
    let up_region = Region::new(&f.chr_up, up_start, up_start + up_len);
    let final_down = &fused[fused.len() - 1];
    let down_len = final_down.gene_end - final_down.gene_start;
    let down_start = if f.strand_down < 0 { f.pos_down - down_len } else { f.pos_down };
    let down_region = Region::new(&f.chr_down, down_start, down_start + down_len);
    for d in domains {
        if d.transcript == f.transcript_up && d.chromosome == up_region.chromosome && up_region.end > d.start && up_region.start < d.end {
            let c = convert_region(f.strand_up, &up_region, d.start, d.end);
            let mut n = d.clone();
            n.chromosome = f.name();
            n.start = i32::max(c.0, first_up.gene_start);
            n.end = i32::min(c.1, first_up.gene_end);
            out.push(n);
        }
        if d.transcript == f.transcript_down && d.chromosome == down_region.chromosome && down_region.end > d.start && down_region.start < d.end {
            let c = convert_region(f.strand_down, &down_region, d.start, d.end);
            let mut n = d.clone();
            n.chromosome = f.name();
            n.start = i32::max(c.0 + first_up.gene_end, final_down.gene_start);
            n.end = i32::min(c.1 + first_up.gene_end, final_down.gene_end);
            out.push(n);
        }
    }
    out
}

/// ScaleIntrons: introns shortened to 10 bases
pub(crate) struct ScaleIntrons {
    pub(crate) introns: Vec<Region>,
}

impl ScaleIntrons {
    fn from_exons(exons: &[FusedExon]) -> ScaleIntrons {
        let regions: Vec<Region> = exons.iter().map(|e| Region::new(&e.chromosome, e.start, e.end)).collect();
        Self::from_regions(&regions)
    }

    /// ScaleIntrons.introns: the gaps between consecutive regions
    pub(crate) fn from_regions(exons: &[Region]) -> ScaleIntrons {
        let mut introns = Vec::new();
        for i in 0..exons.len().saturating_sub(1) {
            let (cur, next) = (&exons[i], &exons[i + 1]);
            if next.start > cur.end {
                introns.push(Region::new(&cur.chromosome, cur.end + 1, next.start - 1));
            }
        }
        ScaleIntrons { introns }
    }

    /// scaleIntronsFromRegion
    pub(crate) fn scale_region(&self, r: &Region) -> Region {
        Region::new(&r.chromosome, r.start + self.offset(r.start as i64), r.end + self.offset(r.end as i64))
    }

    fn offset(&self, position: i64) -> i32 {
        let mut offset: i64 = 0;
        for intron in &self.introns {
            if position > intron.start as i64 {
                let bases = intron.bases() as i64;
                let distance = i64::min(bases, position - intron.start as i64 + 1);
                let proportion = distance as f64 / bases as f64;
                offset += java_round(proportion * 10.0) - distance;
            }
        }
        offset as i32
    }
}

/// java.awt.Color.HSBtoRGB, in float arithmetic as Java does it
pub fn hsb_to_rgb(hue: f32, saturation: f32, brightness: f32) -> (u8, u8, u8) {
    let c = |v: f32| (v * 255.0f32 + 0.5f32) as i32 as u8;
    if saturation == 0.0 {
        let v = c(brightness);
        return (v, v, v);
    }
    let h = (hue - hue.floor()) * 6.0f32;
    let f = h - h.floor();
    let p = brightness * (1.0f32 - saturation);
    let q = brightness * (1.0f32 - saturation * f);
    let t = brightness * (1.0f32 - (saturation * (1.0f32 - f)));
    match h as i32 {
        0 => (c(brightness), c(t), c(p)),
        1 => (c(q), c(brightness), c(p)),
        2 => (c(p), c(brightness), c(t)),
        3 => (c(p), c(q), c(brightness)),
        4 => (c(t), c(p), c(brightness)),
        5 => (c(brightness), c(p), c(q)),
        _ => (0, 0, 0),
    }
}

/// ProteinDomainColors.hue
pub fn protein_hue(i: i32, max_domains: i32) -> f32 {
    let colour_distance = 1f32 / (((max_domains + 2) as f32 / 2f32) as f64).ceil() as f32 / 2f32;
    if i % 2 == 0 {
        21f32 / 360f32 + (i / 2 + 1) as f32 * colour_distance
    } else {
        201f32 / 360f32 + ((i as f32 / 2f32) as f64).ceil() as f32 * colour_distance
    }
}

/// ProteinDomainColors: UTR grey, the other domains spread over hues (in the iteration order of Java's
/// HashSet of domain names)
pub struct ProteinDomainColours {
    map: Vec<(String, (u8, u8, u8))>,
}

impl ProteinDomainColours {
    pub fn new(domains: &[VisProteinDomain]) -> ProteinDomainColours {
        let mut names: JavaHashSet<String> = JavaHashSet::new();
        for d in domains {
            names.insert(d.name().to_string());
        }
        let mut map = vec![(PROTEIN_DOMAIN_UTR.to_string(), (128u8, 128u8, 128u8))];
        let new: Vec<String> = names.iter().filter(|n| n.as_str() != PROTEIN_DOMAIN_UTR).cloned().collect();
        for (i, n) in new.iter().enumerate() {
            map.push((n.clone(), hsb_to_rgb(protein_hue(i as i32, new.len() as i32), 0.5f32, 0.8392f32)));
        }
        ProteinDomainColours { map }
    }

    pub fn colour(&self, d: &VisProteinDomain) -> (u8, u8, u8) {
        if d.name() == PROTEIN_DOMAIN_UTR {
            return if d.start == 1 { (107, 174, 214) } else { (214, 144, 107) };
        }
        self.map.iter().find(|(n, _)| n == d.name()).map(|(_, c)| *c).unwrap_or((0, 0, 0))
    }
}

pub struct FusionData {
    pub exons: Vec<FusedExon>,
    pub domains: Vec<VisProteinDomain>,
    pub colours: ProteinDomainColours,
}

impl FusionData {
    pub fn new(fusions: &[VisFusion], exons: &[VisGeneExon], domains: &[VisProteinDomain]) -> Result<FusionData> {
        let exonic = exonic_protein_domains(domains, exons);
        let mut final_exons = Vec::new();
        let mut final_domains = Vec::new();
        for f in fusions {
            let fe = fused_exons(f, exons);
            let fd = fused_protein_domains(f, &fe, &exonic);
            let scaler = ScaleIntrons::from_exons(&fe);
            let fe: Vec<FusedExon> = fe
                .into_iter()
                .map(|mut e| {
                    let (s, en, gs, ge) = (e.start, e.end, e.gene_start, e.gene_end);
                    e.start = s + scaler.offset(s as i64);
                    e.end = en + scaler.offset(en as i64);
                    e.gene_start = gs + scaler.offset(gs as i64);
                    e.gene_end = ge + scaler.offset(ge as i64);
                    e
                })
                .collect();
            let fd: Vec<VisProteinDomain> = fd
                .into_iter()
                .map(|mut d| {
                    let (s, e) = (d.start, d.end);
                    d.start = s + scaler.offset(s as i64);
                    d.end = e + scaler.offset(e as i64);
                    d
                })
                .collect();
            let mut positions = Vec::new();
            for e in &fe {
                for p in [e.gene_start, e.gene_end, e.start, e.end] {
                    positions.push(GenomePosition::new(&e.fusion, p));
                }
            }
            for d in &fd {
                positions.push(GenomePosition::new(&d.chromosome, d.start));
                positions.push(GenomePosition::new(&d.chromosome, d.end));
            }
            let sp = ScalePosition::new(&positions);
            for mut e in fe {
                let c = sp.contig(&e.fusion)?;
                e.gene_start = c.scale(e.gene_start)?;
                e.gene_end = c.scale(e.gene_end)?;
                e.start = c.scale(e.start)?;
                e.end = c.scale(e.end)?;
                final_exons.push(e);
            }
            for mut d in fd {
                let c = sp.contig(&d.chromosome)?;
                d.start = c.interpolate(d.start);
                d.end = c.interpolate(d.end);
                final_domains.push(d);
            }
        }
        let colours = ProteinDomainColours::new(&final_domains);
        Ok(FusionData { exons: final_exons, domains: final_domains, colours })
    }

    /// FusedExons.write's lines
    pub fn exon_lines(&self) -> Vec<String> {
        let mut out = vec!["sampleId\tclusterId\tfusion\tgene\tgeneStart\tgeneEnd\tisUpstream\tchromosome\tstart\tend\trank\tskipped\ttranscript".to_string()];
        for e in &self.exons {
            out.push(
                [
                    String::new(),
                    e.cluster_id.to_string(),
                    e.fusion.clone(),
                    e.gene.clone(),
                    e.gene_start.to_string(),
                    e.gene_end.to_string(),
                    e.is_upstream.to_string(),
                    e.chromosome.clone(),
                    e.start.to_string(),
                    e.end.to_string(),
                    e.rank.to_string(),
                    e.skipped.to_string(),
                    e.transcript.clone(),
                ]
                .join("\t"),
            );
        }
        out
    }

    /// FusedProteinDomains.write's lines
    pub fn domain_lines(&self) -> Vec<String> {
        let mut out = vec!["sampleId\tclusterId\tfusion\tstart\tend\tname\tcolor\ttranscript".to_string()];
        for d in &self.domains {
            let c = self.colours.colour(d);
            out.push(
                [
                    String::new(),
                    d.cluster_id.to_string(),
                    d.chromosome.clone(),
                    d.start.to_string(),
                    d.end.to_string(),
                    d.name().to_string(),
                    format!("#{:02X}{:02X}{:02X}", c.0, c.1, c.2),
                    d.transcript.clone(),
                ]
                .join("\t"),
            );
        }
        out
    }
}
