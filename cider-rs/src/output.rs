//! The output files: <sample>.cider.vdj.tsv.gz (VDJSequenceTsvWriter), .cider.alignment_match.tsv.gz
//! (AlignmentMatchTsvWriter), .cider.locus_stats.tsv (CiderLocusStatsWriter), .cider.layout.gz
//! (VJReadLayoutFile) and .cider.primer_match.tsv (VdjPrimerMatchTsv). Gzip as java.util.zip writes it
//! (zlib level 6), so the compressed files are Java's bytes too.

use std::collections::HashSet;
use std::io::Write;

use anyhow::{Context, Result};
use hmf_common::codon::{amino_acid_from_bases, STOP_AMINO_ACID};
use hmf_common::csv::CsvPrinter;
use hmf_common::fmt::java_double_to_string;

use crate::alignment::{Alignment, AlignmentAnnotation, AlignmentStatus, GeneAnnotation};
use crate::annotate::{locus, split_read_counts, VdjAnnotation};
use crate::genes::{IgTcrGene, IgTcrLocus, VjGeneType};
use crate::layout::{counts_to_string, ReadLayout};
use crate::layout_builder::{candidate, LayoutBuildResult, VjReadLayoutBuilder};
use crate::primer::VdjPrimerMatch;
use crate::read::cigar_to_string;
use crate::vdj::{insert_dashes, VdjSequence};

/// a plain or gzipped file (FileWriterUtils.createBufferedWriter / createGzipBufferedWriter)
pub enum Sink {
    Plain(std::fs::File),
    Gz(flate2::write::GzEncoder<std::fs::File>),
}

impl Write for Sink {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Sink::Plain(f) => f.write(buf),
            Sink::Gz(g) => g.write(buf),
        }
    }
    /// a flush of the gzip stream would add an empty block (Z_SYNC_FLUSH), which Java does not write
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Sink::Plain(f) => f.flush(),
            Sink::Gz(_) => Ok(()),
        }
    }
}

pub type TextWriter = std::io::BufWriter<Sink>;

pub fn create_writer(path: &str, gzip: bool) -> Result<TextWriter> {
    let file = std::fs::File::create(path).with_context(|| format!("creating {path}"))?;
    let sink = if gzip { Sink::Gz(flate2::write::GzEncoder::new(file, flate2::Compression::new(6))) } else { Sink::Plain(file) };
    Ok(std::io::BufWriter::with_capacity(1 << 16, sink))
}

pub fn finish(w: TextWriter, path: &str) -> Result<()> {
    let sink = w.into_inner().map_err(|e| anyhow::anyhow!("writing {path}: {}", e.error()))?;
    match sink {
        Sink::Plain(mut f) => f.flush().with_context(|| format!("writing {path}"))?,
        Sink::Gz(g) => {
            g.finish().with_context(|| format!("writing {path}"))?;
        }
    }
    Ok(())
}

/// CiderFormatter.aminoAcidFromBases: stops as '_'
pub fn format_aa(dna: &[u8]) -> String {
    amino_acid_from_bases(dna).replace(STOP_AMINO_ACID as char, "_")
}

/// CiderFormatter.cdr3AminoAcid: "fs" added when out of frame
pub fn cdr3_amino_acid(vdj: &VdjSequence) -> String {
    format!("{}{}", format_aa(vdj.cdr3_sequence()), if vdj.is_in_frame() { "" } else { "fs" })
}

fn s(b: &[u8]) -> String {
    String::from_utf8_lossy(b).into_owned()
}

const VDJ_COLUMNS: &[&str] = &[
    "cdr3Seq",
    "cdr3AA",
    "locus",
    "filter",
    "alignmentStatus",
    "minHighQualBaseReads",
    "assignedReads",
    "vAlignedReads",
    "jAlignedReads",
    "inFrame",
    "containsStop",
    "vType",
    "vAnchorStart",
    "vAnchorEnd",
    "vAnchorSeq",
    "vAnchorTemplateSeq",
    "vAnchorAA",
    "vAnchorTemplateAA",
    "vMatchMethod",
    "vSimilarityScore",
    "vNonSplitReads",
    "jType",
    "jAnchorStart",
    "jAnchorEnd",
    "jAnchorSeq",
    "jAnchorTemplateSeq",
    "jAnchorAA",
    "jAnchorTemplateAA",
    "jMatchMethod",
    "jSimilarityScore",
    "jNonSplitReads",
    "vGene",
    "vGeneSupplementary",
    "vPIdent",
    "vPIdentLen",
    "vPIdentIndelBases",
    "vPIdentClipBases",
    "vImgtRefUnaligned",
    "SHMStatus",
    "vAlignStart",
    "vAlignEnd",
    "dGene",
    "dGeneSupplementary",
    "dAlignStart",
    "dAlignEnd",
    "jGene",
    "jGeneSupplementary",
    "jPIdent",
    "jPIdentLen",
    "jPIdentIndelBases",
    "jPIdentClipBases",
    "jImgtRefUnaligned",
    "jAlignStart",
    "jAlignEnd",
    "vPrimerMatches",
    "jPrimerMatches",
    "layoutId",
    "fullSeq",
    "support",
];

/// formatGeneSupplementary: the other alleles, distinct and sorted, joined by ';'
fn gene_supplementary(g: Option<&GeneAnnotation>) -> Option<String> {
    let g = g?;
    let primary = g.gene.gene_allele();
    let mut set: Vec<String> = Vec::new();
    for s in g.supplementary_genes.iter().map(IgTcrGene::gene_allele).filter(|a| *a != primary) {
        if !set.contains(&s) {
            set.push(s);
        }
    }
    set.sort();
    Some(set.join(";"))
}

fn opt<T: ToString>(v: Option<T>) -> Option<String> {
    v.map(|x| x.to_string())
}

/// VDJSequenceTsvWriter.writeVDJSequences
pub fn write_vdj_sequences(path: &str, vdjs: &[VdjSequence], annotations: &[VdjAnnotation], alignments: &[AlignmentAnnotation]) -> Result<()> {
    let w = create_writer(path, true)?;
    let mut p = CsvPrinter::new(w, b'\t');
    p.print_record(VDJ_COLUMNS)?;
    for a in annotations {
        let vdj = &vdjs[a.vdj];
        let al = a.alignment.map(|k| &alignments[k]);
        let offset = al.map(|x| x.query_range.0).unwrap_or(0);
        let v = vdj.v_anchor.as_ref();
        let j = vdj.j_anchor.as_ref();
        let vg = al.and_then(|x| x.v_gene.as_ref());
        let dg = al.and_then(|x| x.d_gene.as_ref());
        let jg = al.and_then(|x| x.j_gene.as_ref());
        let vc = vg.and_then(|g| g.comparison.as_ref());
        let jc = jg.and_then(|g| g.comparison.as_ref());
        let align_start = |g: Option<&GeneAnnotation>| g.map(|g| (g.alignment.query_start + offset).to_string());
        let align_end = |g: Option<&GeneAnnotation>| g.map(|g| (g.alignment.query_end + offset + 1).to_string());
        let values: Vec<Option<String>> = vec![
            Some(s(vdj.cdr3_sequence())),
            Some(cdr3_amino_acid(vdj)),
            Some(locus(vdj).pretty_print().to_string()),
            Some(a.filters.iter().map(|f| f.name()).collect::<Vec<_>>().join(";")),
            Some(al.map(|x| x.status).unwrap_or(AlignmentStatus::SkippedAlign).name().to_string()),
            Some(a.cdr3_support_min.to_string()),
            Some(vdj.num_reads().to_string()),
            Some(a.v_aligned_reads.to_string()),
            Some(a.j_aligned_reads.to_string()),
            Some(vdj.is_in_frame().to_string()),
            Some(vdj.amino_acid_sequence().contains(STOP_AMINO_ACID as char).to_string()),
            v.map(|x| x.gene_type.name().to_string()),
            v.map(|_| vdj.layout_slice_start.to_string()),
            v.map(|x| (vdj.layout_slice_start + x.anchor_boundary).to_string()),
            vdj.v_anchor_sequence().map(s),
            v.map(|x| x.template_anchor_seq.clone()),
            vdj.v_anchor_sequence().map(|seq| format_aa(&seq[seq.len() % 3..])),
            v.map(|x| format_aa(x.template_anchor_seq.as_bytes())),
            v.map(|x| x.match_method().to_string()),
            opt(a.v_similarity_score),
            Some(a.v_non_split_reads.to_string()),
            j.map(|x| x.gene_type.name().to_string()),
            j.map(|x| (vdj.layout_slice_start + x.anchor_boundary).to_string()),
            j.map(|_| vdj.layout_slice_end.to_string()),
            vdj.j_anchor_sequence().map(s),
            j.map(|x| x.template_anchor_seq.clone()),
            vdj.j_anchor_sequence().map(format_aa),
            j.map(|x| format_aa(x.template_anchor_seq.as_bytes())),
            j.map(|x| x.match_method().to_string()),
            opt(a.j_similarity_score),
            Some(a.j_non_split_reads.to_string()),
            vg.map(|g| g.gene.gene_allele()),
            gene_supplementary(vg),
            vc.map(|c| java_double_to_string(c.pct_identity)),
            opt(vc.map(|c| c.seq_length)),
            opt(vc.map(|c| c.indel_bases)),
            opt(vc.map(|c| c.clip_bases)),
            opt(vc.map(|c| c.imgt_ref_unaligned)),
            a.somatic_hypermutation_status.map(|x| x.to_string()),
            align_start(vg),
            align_end(vg),
            dg.map(|g| g.gene.gene_allele()),
            gene_supplementary(dg),
            align_start(dg),
            align_end(dg),
            jg.map(|g| g.gene.gene_allele()),
            gene_supplementary(jg),
            jc.map(|c| java_double_to_string(c.pct_identity)),
            opt(jc.map(|c| c.seq_length)),
            opt(jc.map(|c| c.indel_bases)),
            opt(jc.map(|c| c.clip_bases)),
            opt(jc.map(|c| c.imgt_ref_unaligned)),
            align_start(jg),
            align_end(jg),
            Some(a.v_primer_match_count.to_string()),
            Some(a.j_primer_match_count.to_string()),
            Some(vdj.layout.id.clone()),
            Some(vdj.layout.consensus_string()),
            Some(counts_to_string(&vdj.layout.high_qual_support_counts())),
        ];
        for v in values {
            p.print_opt(v)?;
        }
        p.println()?;
    }
    finish(p.into_inner(), path)
}

/// AlignmentMatchTsvWriter.write: a row per V, D, J and full alignment, "null" for no value
pub fn write_alignment_matches(path: &str, vdjs: &[VdjSequence], annotations: &[AlignmentAnnotation]) -> Result<()> {
    let w = create_writer(path, true)?;
    let mut p = CsvPrinter::new(w, b'\t').with_null_string("null");
    p.print_record(&[
        "cdr3Seq",
        "fullSeq",
        "matchType",
        "gene",
        "functionality",
        "layoutAlignStart",
        "layoutAlignEnd",
        "alignScore",
        "refContig",
        "refStart",
        "refEnd",
        "strand",
        "refContigLength",
        "cigar",
        "editDistance",
        "querySeqStart",
        "querySeqEnd",
        "querySeq",
    ])?;
    for a in annotations {
        let vdj = &vdjs[a.vdj];
        let rows: [(&str, Option<&IgTcrGene>, Option<&Alignment>); 4] = [
            ("V", a.v_gene.as_ref().map(|g| &g.gene), a.v_gene.as_ref().map(|g| &g.alignment)),
            ("D", a.d_gene.as_ref().map(|g| &g.gene), a.d_gene.as_ref().map(|g| &g.alignment)),
            ("J", a.j_gene.as_ref().map(|g| &g.gene), a.j_gene.as_ref().map(|g| &g.alignment)),
            ("full", None, a.full_alignment.as_ref()),
        ];
        for (kind, gene, alignment) in rows {
            let Some(al) = alignment else { continue };
            let offset = a.query_range.0;
            let values: Vec<Option<String>> = vec![
                Some(s(vdj.cdr3_sequence())),
                Some(vdj.layout.consensus_string()),
                Some(kind.to_string()),
                gene.map(|g| g.gene_allele()),
                gene.map(|g| g.functionality.to_code().to_string()),
                Some((offset + al.query_start).to_string()),
                Some((offset + al.query_end + 1).to_string()),
                Some(al.alignment_score.to_string()),
                Some(al.ref_contig.clone()),
                Some((al.ref_start + 1).to_string()),
                Some((al.ref_end + 1).to_string()),
                Some(al.strand.as_char().to_string()),
                Some(al.ref_contig_length.to_string()),
                Some(if al.cigar.is_empty() { String::new() } else { cigar_to_string(&al.cigar) }),
                Some(al.edit_distance.to_string()),
                Some(a.query_range.0.to_string()),
                Some((a.query_range.1 + 1).to_string()),
                Some(al.query_seq.clone()),
            ];
            for v in values {
                p.print_opt(v)?;
            }
            p.println()?;
        }
    }
    finish(p.into_inner(), path)
}

/// CiderLocusStatsWriter.writeLocusStats
pub fn write_locus_stats(
    path: &str,
    results: &[(VjGeneType, LayoutBuildResult)],
    vdjs: &[VdjSequence],
    annotations: &[VdjAnnotation],
) -> Result<()> {
    let w = create_writer(path, false)?;
    let mut p = CsvPrinter::new(w, b'\t');
    p.print_record(&["locus", "readsUsed", "readsTotal", "downSampled", "sequences", "passSequences"])?;
    for l in IgTcrLocus::ALL {
        let of_locus: Vec<&LayoutBuildResult> = results.iter().filter(|(t, _)| t.locus() == l).map(|(_, r)| r).collect();
        let total: usize = of_locus.iter().map(|r| r.num_candidate_reads).sum();
        let mut keys = HashSet::new();
        for r in &of_locus {
            for layout in &r.layouts {
                for read in layout.reads() {
                    keys.insert(read.key.clone());
                }
            }
        }
        let down_sampled = of_locus.iter().any(|r| r.down_sampled);
        let of = annotations.iter().filter(|a| locus(&vdjs[a.vdj]) == l);
        let (n, pass) = of.fold((0, 0), |(n, p), a| (n + 1, p + a.passes_filter() as i32));
        for v in [l.pretty_print().to_string(), keys.len().to_string(), total.to_string(), down_sampled.to_string(), n.to_string(), pass.to_string()] {
            p.print(&v)?;
        }
        p.println()?;
    }
    finish(p.into_inner(), path)
}

const ANSI_RESET: &str = "\u{1b}[0m";

/// IGV colours: A green, C purple, G yellow, T red
fn colorise(s: &str) -> String {
    let mut out = String::with_capacity(s.len() * 6);
    for c in s.chars() {
        match c {
            'A' => out.push_str("\u{1b}[32m"),
            'C' => out.push_str("\u{1b}[35m"),
            'G' => out.push_str("\u{1b}[33m"),
            'T' => out.push_str("\u{1b}[31m"),
            _ => {}
        }
        out.push(c);
    }
    out.push_str(ANSI_RESET);
    out
}

/// VJReadLayoutFile.writeLayouts: by gene type, each type's layouts by descending read count
pub fn write_layouts(path: &str, results: &[(VjGeneType, LayoutBuildResult)]) -> Result<()> {
    let mut w = create_writer(path, true)?;
    let mut sorted: Vec<&(VjGeneType, LayoutBuildResult)> = results.iter().collect();
    sorted.sort_by_key(|(t, _)| *t);
    for (gene_type, r) in sorted {
        let mut layouts: Vec<&ReadLayout> = r.layouts.iter().map(|l| l.as_ref()).collect();
        layouts.sort_by(|a, b| b.reads().len().cmp(&a.reads().len()));
        for layout in layouts {
            write_layout(&mut w, layout, *gene_type)?;
        }
    }
    finish(w, path)
}

fn write_layout(w: &mut TextWriter, layout: &ReadLayout, gene_type: VjGeneType) -> Result<()> {
    let (split5, split10) = split_read_counts(layout);
    let (first, last) = VjReadLayoutBuilder::extrapolated_anchor_range(gene_type.vj(), layout);
    let sequence = layout.consensus_string();
    let support = layout.high_qual_support_string();
    let shift = first.rem_euclid(3) as usize;
    let codon_aligned = &sequence.as_bytes()[shift.min(sequence.len())..];
    let aa = insert_dashes(&amino_acid_from_bases(codon_aligned), &[first / 3, (last + 1) / 3 - 1 + 1]);
    let dashes = |x: &str| insert_dashes(x, &[first, last + 1]);
    writeln!(
        w,
        "{} type: {}, read count: {}, split(5) read count: {}, split(10) read count: {}, AA: {}",
        layout.id,
        gene_type.name(),
        layout.reads().len(),
        split5,
        split10,
        aa
    )?;
    writeln!(w, "    {}", colorise(&dashes(&sequence)))?;
    writeln!(w, "    {}", dashes(&support))?;
    for r in layout.reads() {
        let c = candidate(r);
        let pad = (layout.aligned_position - r.aligned_position).max(0) as usize;
        let seq = format!("{}{}", " ".repeat(pad), String::from_utf8_lossy(&r.sequence));
        let qual: String = format!("{}{}", " ".repeat(pad), r.base_qualities.iter().map(|&q| (q + 33) as char).collect::<String>());
        writeln!(w, "    read: {}", c.read.to_hts_string())?;
        writeln!(w, "    {}", colorise(&dashes(&seq)))?;
        writeln!(w, "    {}", dashes(&qual))?;
    }
    Ok(())
}

/// VdjPrimerMatchTsv.writePrimerMatches
pub fn write_primer_matches(path: &str, vdjs: &[VdjSequence], matches: &[VdjPrimerMatch]) -> Result<()> {
    let w = create_writer(path, false)?;
    let mut p = CsvPrinter::new(w, b'\t');
    p.print_record(&["cdr3Seq", "cdr3AA", "primerGene", "primerName", "numMismatch", "primerSequence", "vdjFullSeq"])?;
    for m in matches {
        let vdj = &vdjs[m.vdj];
        for v in [
            s(vdj.cdr3_sequence()),
            cdr3_amino_acid(vdj),
            m.primer.target.clone(),
            m.primer.name.clone(),
            m.num_mismatch.to_string(),
            m.primer.sequence.clone(),
            m.full_vdj_sequence.clone(),
        ] {
            p.print(&v)?;
        }
        p.println()?;
    }
    finish(p.into_inner(), path)
}
