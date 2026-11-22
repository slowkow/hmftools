//! VcfWriter: the raw breakend VCF (bgzip compressed with a tabix index), laid out as htsjdk writes it.

use std::io::Write;

use anyhow::{Context, Result};
use rust_htslib::bgzf;

use crate::assembly::alignment::alternative_alignment::AlternativeAlignment;
use crate::assembly::alignment::breakend::{BreakendP, BreakendSupport};
use crate::assembly::common_utils::is_valid_dna_base;
use crate::assembly::config::{AssemblyConfig, WriteType};
use crate::assembly::constants::BREAKEND_REQ_VALID_FRAGMENT_LENGTH_PERC;
use crate::assembly::line_utils::is_mobile_line_element;
use crate::assembly::orientation::Orientation;
use crate::common::fmt::vcf_double;

/// The version of the Java tool this is a port of, as written to the header.
const ESVEE_VERSION: &str = "2.1";

/// Header lines in the order htsjdk sorts them, between the file format line and the contig lines.
const HEADER_LINES: &str = r#"##FILTER=<ID=PON,Description="Found in panel of normals">
##FILTER=<ID=dedup,Description="Event is duplicate of another">
##FILTER=<ID=delShortLowVaf,Description="Short low-VAF DEL artefacts with homology">
##FILTER=<ID=invShortFragLowVaf,Description="Short (<300b) low-VAF INV artefacts">
##FILTER=<ID=invShortIsolated,Description="Short (<100b) isolated INV artefacts">
##FILTER=<ID=invShortLowVafHom,Description="Short (<3K) low-VAF INV artefacts with homology">
##FILTER=<ID=lineSource,Description="Remote LINE source">
##FILTER=<ID=minAF,Description="Variant allele fraction too low">
##FILTER=<ID=minAnchorLength,Description="Trimmed anchor length is too short">
##FILTER=<ID=minLength,Description="Variant is too short">
##FILTER=<ID=minQual,Description="Insufficient quality">
##FILTER=<ID=minSupport,Description="Below minimum fragment support">
##FILTER=<ID=offTarget,Description="Outside target panel">
##FILTER=<ID=sbxArtefact,Description="SBX artefact">
##FILTER=<ID=sbxInvZeroMotif,Description="SBX artefact, zero-length INV with known motif">
##FILTER=<ID=sbxStrandBias,Description="SBX strand-bias">
##FILTER=<ID=sgl,Description="SGLs filtered entirely">
##FILTER=<ID=sglStrandBias,Description="Single breakend with excessive strand bias and not a mobile line insertion">
##FILTER=<ID=shortFrags,Description="Average variant fragment length is too short">
##FILTER=<ID=strandBias,Description="Outside valid strand bias range">
##FILTER=<ID=unpairedPositionRange,Description="Unpaired position range">
##FILTER=<ID=weakJuncMinAF,Description="Low allele fraction at weak junction">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DF,Number=1,Type=Integer,Description="Count of discordant fragments with a read either side of the breakend">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=QUAL,Number=1,Type=Integer,Description="Variant quality">
##FORMAT=<ID=SB,Number=1,Type=Float,Description="Strand read bias">
##FORMAT=<ID=SF,Number=1,Type=Integer,Description="Count of fragments supporting the breakend with a read overlapping the breakend">
##FORMAT=<ID=VF,Number=1,Type=Integer,Description="Total variant fragments supporting the breakend">
##INFO=<ID=ALTALN,Number=1,Type=String,Description="Alternative alignments for low map-qual breakends">
##INFO=<ID=ASMID,Number=.,Type=Integer,Description="Unique id(s) of assembly(s) containing the breakend">
##INFO=<ID=ASMINFO,Number=1,Type=String,Description="Original junction coordinates">
##INFO=<ID=ASMLEN,Number=1,Type=Integer,Description="Assembly sequence length">
##INFO=<ID=ASMLNKS,Number=.,Type=Integer,Description="Id(s) of breakend(s) linked by assembly">
##INFO=<ID=ASMSEG,Number=.,Type=Integer,Description="Segment indices in assembly(s) containing the breakend">
##INFO=<ID=AVGLEN,Number=1,Type=Integer,Description="Average variant fragment length">
##INFO=<ID=BEAOR,Number=.,Type=Integer,Description="Breakend orientation(s) in assembly(s)">
##INFO=<ID=BEAPOS,Number=.,Type=Integer,Description="Breakend position(s) in assembly(s)">
##INFO=<ID=BEOR,Number=.,Type=Integer,Description="Breakend orientation(s) in reference genome">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around position">
##INFO=<ID=DF,Number=1,Type=Integer,Description="Count of discordant fragments with a read either side of the breakend">
##INFO=<ID=HOMSEQ,Number=1,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=IHOMPOS,Number=2,Type=Integer,Description="Position of inexact homology">
##INFO=<ID=INSALN,Number=1,Type=String,Description="Alternative alignments of insert sequence">
##INFO=<ID=LINE,Number=0,Type=Flag,Description="LINE insertion site">
##INFO=<ID=MATEID,Number=1,Type=String,Description="Mate breakend ID">
##INFO=<ID=MLR,Number=1,Type=Integer,Description="Max local indel repeat around breakend">
##INFO=<ID=PJRR,Number=1,Type=Float,Description="Proximate junction read ratio">
##INFO=<ID=SAGA,Number=1,Type=String,Description="ID and breakends of the matched SAGA variant">
##INFO=<ID=SAGAINF,Number=0,Type=Flag,Description="Indicates the breakend is inferred from the matched SAGA variant">
##INFO=<ID=SEGALEN,Number=1,Type=Integer,Description="Aligned length of segment(s) in reference genome">
##INFO=<ID=SEGID,Number=.,Type=Integer,Description="Unique id(s) of segment(s) containing the breakend">
##INFO=<ID=SEGMAPQ,Number=1,Type=Integer,Description="MAPQ of segment containing the breakend with highest QUAL contribution">
##INFO=<ID=SEGRL,Number=1,Type=Integer,Description="Repeat length of segment with highest QUAL contribution">
##INFO=<ID=SEGSCO,Number=1,Type=Integer,Description="Alignment score of segments containing the breakend with highest QUAL contribution">
##INFO=<ID=SF,Number=1,Type=Integer,Description="Count of fragments supporting the breakend with a read overlapping the breakend">
##INFO=<ID=SRRR,Number=1,Type=Float,Description="Supplementary remote region ratio">
##INFO=<ID=SVID,Number=1,Type=String,Description="ID shared by two breakends">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=TPR,Number=2,Type=Integer,Description="Unpaired read three-prime position range">
##INFO=<ID=UFP,Number=1,Type=Integer,Description="Distinct fragment positions">
##INFO=<ID=VF,Number=1,Type=Integer,Description="Total variant fragments supporting the breakend">
"#;

pub struct VcfWriter<'a> {
    config: &'a AssemblyConfig,
    path: String,
    out: Option<bgzf::Writer>,
    /// the index into each breakend's sample support for each written sample (reference samples first)
    sample_support_indices: Vec<usize>,
    variants: Vec<String>,
}

/// hmf-common VariantAltInsertCoords.formPairedAltString
pub fn form_paired_alt_string(alt: &str, insert_sequence: &str, chromosome: &str, position: i32, orient_start: Orientation, orient_end: Orientation) -> String {
    if orient_start.is_forward() && orient_end.is_reverse() {
        format!("{alt}{insert_sequence}[{chromosome}:{position}[")
    } else if orient_start.is_forward() && orient_end.is_forward() {
        format!("{alt}{insert_sequence}]{chromosome}:{position}]")
    } else if orient_start.is_reverse() && orient_end.is_reverse() {
        format!("[{chromosome}:{position}[{insert_sequence}{alt}")
    } else {
        format!("]{chromosome}:{position}]{insert_sequence}{alt}")
    }
}

/// hmf-common VariantAltInsertCoords.formSingleAltString
pub fn form_single_alt_string(alt: &str, insert_sequence: &str, orientation: Orientation) -> String {
    if orientation.is_forward() { format!("{alt}{insert_sequence}.") } else { format!(".{insert_sequence}{alt}") }
}

/// The genotype column of one sample: GT:AD:DF:DP:SB:SF:VF
pub fn build_genotype(breakend_support: &BreakendSupport) -> String {
    let alt_support = breakend_support.total_support();
    format!(
        "./.:0,{}:{}:0:{}:{}:{}",
        alt_support,
        breakend_support.discordant_fragments,
        vcf_double(breakend_support.strand_bias()),
        breakend_support.split_fragments,
        alt_support
    )
}

fn join_segments<T: ToString>(values: impl Iterator<Item = T>) -> String {
    values.map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}

impl<'a> VcfWriter<'a> {
    pub fn new(config: &'a AssemblyConfig) -> Result<VcfWriter<'a>> {
        let sample_names: Vec<String> = config.reference_ids.iter().chain(config.tumor_ids.iter()).cloned().collect();
        let config_sample_ids = config.combined_sample_ids();
        // a name occurring twice maps to its last index, as with the Java map
        let sample_support_indices: Vec<usize> = sample_names.iter().map(|s| config_sample_ids.iter().rposition(|x| x == s).expect("sample")).collect();
        let mut w = VcfWriter { config, path: String::new(), out: None, sample_support_indices, variants: Vec::new() };
        if config.write_types.contains(&WriteType::Vcf) {
            let path = config.output_filename(WriteType::Vcf);
            let mut out = bgzf::Writer::from_path(&path).with_context(|| format!("creating {path}"))?;
            let mut header = String::new();
            header.push_str("##fileformat=VCFv4.2\n");
            header.push_str(&format!("##EsveeVersion={ESVEE_VERSION}\n"));
            header.push_str(HEADER_LINES);
            if let Some(ref_genome) = config.ref_genome() {
                for (name, length) in ref_genome.chromosomes() {
                    header.push_str(&format!("##contig=<ID={name},length={length}>\n"));
                }
            }
            header.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
            for s in &sample_names {
                header.push('\t');
                header.push_str(s);
            }
            header.push('\n');
            out.write_all(header.as_bytes())?;
            w.path = path;
            w.out = Some(out);
        }
        Ok(w)
    }

    pub fn close(mut self) -> Result<()> {
        let Some(mut out) = self.out.take() else { return Ok(()) };
        for line in &self.variants {
            out.write_all(line.as_bytes())?;
            out.write_all(b"\n")?;
        }
        out.flush()?;
        drop(out);
        let c_path = std::ffi::CString::new(self.path.as_str())?;
        let ret = unsafe { rust_htslib::htslib::tbx_index_build(c_path.as_ptr(), 0, &rust_htslib::htslib::tbx_conf_vcf) };
        if ret != 0 {
            anyhow::bail!("failed to build tabix index for {}", self.path);
        }
        Ok(())
    }

    pub fn add_breakend(&mut self, breakend: BreakendP) {
        let mut genotypes: Vec<String> = Vec::with_capacity(self.sample_support_indices.len());
        let mut total_split_frags = 0;
        let mut total_disc_frags = 0;
        for sample_support_index in &self.sample_support_indices {
            let breakend_support = &breakend.sample_support()[*sample_support_index];
            genotypes.push(build_genotype(breakend_support));
            total_split_frags += breakend_support.split_fragments;
            total_disc_frags += breakend_support.discordant_fragments;
        }
        let qual = breakend.calc_sv_qual();
        let (ref_base, alt_bases) = self.build_allele_info(breakend);
        // INFO attributes in key order
        let mut info: Vec<String> = Vec::with_capacity(32);
        let low_qual_alt_alignments = breakend.low_qual_alt_alignments();
        if !low_qual_alt_alignments.is_empty() {
            info.push(format!("ALTALN={}", AlternativeAlignment::to_vcf_tag(&low_qual_alt_alignments)));
        }
        let assembly_alignment = breakend.assembly();
        let assemblies = assembly_alignment.assemblies();
        let segments = breakend.segments();
        info.push(format!("ASMID={}", assembly_alignment.id()));
        info.push(format!("ASMINFO={}", assembly_alignment.info_limit(10)));
        info.push(format!("ASMLEN={}", assembly_alignment.full_sequence_length()));
        if !breakend.facing_breakends().is_empty() {
            info.push(format!("ASMLNKS={}", join_segments(breakend.facing_breakends().iter().map(|x| x.id()))));
        }
        info.push(format!("ASMSEG={}", join_segments(segments.iter().map(|x| x.index))));
        if breakend.is_single() && breakend.valid_fragment_length_percent() < BREAKEND_REQ_VALID_FRAGMENT_LENGTH_PERC {
            info.push("AVGLEN=0".to_string());
        } else {
            info.push(format!("AVGLEN={}", breakend.average_fragment_length()));
        }
        info.push(format!("BEAOR={}", join_segments(segments.iter().map(|x| x.alignment.orientation().as_byte()))));
        info.push(format!("BEAPOS={}", join_segments(segments.iter().map(|x| x.sequence_index))));
        info.push(format!("BEOR={}", breakend.orient.as_byte()));
        if breakend.homology.exists() {
            info.push(format!("CIPOS={},{}", breakend.homology.exact_start, breakend.homology.exact_end));
        } else {
            info.push("CIPOS=0,0".to_string());
        }
        info.push(format!("DF={total_disc_frags}"));
        if breakend.homology.exists() {
            if !breakend.homology.homology.is_empty() {
                info.push(format!("HOMSEQ={}", breakend.homology.homology));
            }
            info.push(format!("IHOMPOS={},{}", breakend.homology.inexact_start, breakend.homology.inexact_end));
        } else {
            info.push("IHOMPOS=0,0".to_string());
        }
        let alt_alignments = breakend.alternative_alignments();
        if !alt_alignments.is_empty() {
            info.push(format!("INSALN={}", AlternativeAlignment::to_vcf_tag(alt_alignments)));
        }
        if assemblies.iter().any(|x| x.has_line_sequence()) && is_mobile_line_element(breakend.orient, &breakend.inserted_bases) {
            info.push("LINE".to_string());
        }
        if let Some(other) = breakend.other_breakend() {
            info.push(format!("MATEID={}", other.id()));
        }
        if breakend.max_local_repeat() > 0 {
            info.push(format!("MLR={}", breakend.max_local_repeat()));
        }
        let mut min_proximate_junc_ratio = -1.0;
        let mut min_supp_remote_region_ratio = -1.0;
        for assembly in &assemblies {
            if assembly.junction().soft_clip_based() {
                if min_proximate_junc_ratio < 0.0 || assembly.stats().proximate_junc_read_ratio < min_proximate_junc_ratio {
                    min_proximate_junc_ratio = assembly.stats().proximate_junc_read_ratio;
                }
                let supp_remote_region_ratio = assembly.stats().supp_remote_region_ratio();
                if min_supp_remote_region_ratio < 0.0 || supp_remote_region_ratio < min_supp_remote_region_ratio {
                    min_supp_remote_region_ratio = supp_remote_region_ratio;
                }
            }
        }
        if min_proximate_junc_ratio >= 0.0 {
            info.push(format!("PJRR={}", vcf_double(min_proximate_junc_ratio)));
        }
        if let Some(saga_match) = assembly_alignment.saga_match() {
            info.push(format!("SAGA={}", saga_match.variant_no_space_string()));
        }
        if breakend.is_saga_inferred() {
            info.push("SAGAINF".to_string());
        }
        info.push(format!("SEGALEN={}", join_segments(segments.iter().map(|x| x.alignment.aligned_bases()))));
        info.push(format!("SEGID={}", join_segments(segments.iter().map(|x| x.unique_id()))));
        info.push(format!("SEGMAPQ={}", segments.iter().map(|x| x.alignment.map_qual()).max().unwrap_or(0)));
        info.push(format!("SEGRL={}", segments.iter().map(|x| x.alignment.adjusted_alignment()).max().unwrap_or(0)));
        info.push(format!("SEGSCO={}", segments.iter().map(|x| x.alignment.score()).max().unwrap_or(0)));
        info.push(format!("SF={total_split_frags}"));
        if min_supp_remote_region_ratio >= 0.0 {
            info.push(format!("SRRR={}", vcf_double(min_supp_remote_region_ratio)));
        }
        if let Some(other) = breakend.other_breakend() {
            let (lower, upper) = if other.id() > breakend.id() { (breakend.id(), other.id()) } else { (other.id(), breakend.id()) };
            info.push(format!("SVID={lower}_{upper}"));
        }
        info.push(format!("SVTYPE={}", breakend.sv_type().name()));
        let read_orientation_range = breakend.read_orientation_range();
        info.push(format!("TPR={},{}", read_orientation_range[0], read_orientation_range[1]));
        info.push(format!("UFP={}", breakend.unique_fragment_position_count()));
        info.push(format!("VF={}", total_split_frags + total_disc_frags));

        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t.\t{}\tGT:AD:DF:DP:SB:SF:VF\t{}",
            breakend.chromosome,
            breakend.position,
            breakend.id(),
            ref_base as char,
            alt_bases,
            qual,
            info.join(";"),
            genotypes.join("\t")
        );
        self.variants.push(line);
    }

    fn build_allele_info(&self, breakend: BreakendP) -> (u8, String) {
        // a null pointer exception in Java without a reference genome or for a position outside it
        let ref_bases = self.config.ref_genome().expect("ref genome").get_bases(&breakend.chromosome, breakend.position, breakend.position).expect("ref base");
        let mut ref_base = ref_bases[0];
        if !is_valid_dna_base(ref_base) {
            ref_base = b'N';
        }
        let alt_base = (ref_base as char).to_string();
        let alt_bases = match breakend.other_breakend() {
            None => form_single_alt_string(&alt_base, &breakend.inserted_bases, breakend.orient),
            Some(other) => form_paired_alt_string(&alt_base, &breakend.inserted_bases, &other.chromosome, other.position, breakend.orient, other.orient),
        };
        (ref_base, alt_bases)
    }
}
