//! NeoEpitopeWriter: candidate fusions for neo-epitope prediction, written as they are found to
//! <sample>.linx.neoepitope.tsv.

use std::io::Write;

use anyhow::{Context, Result};
use hmf_common::ensembl::BIOTYPE_NONSENSE_MED_DECAY;
use hmf_common::fmt::java_fixed;
use hmf_common::known_fusion::KnownFusionType::*;
use hmf_common::known_fusion::{switch_stream, KnownFusionCache, FS_DOWN, FS_UP, POS_STRAND};

use super::disruption::DisruptionFinder;
use super::*;
use crate::gene::{BreakendGene, BreakendTrans, GeneCache, TranscriptRegionType};
use crate::types::*;

struct NeoEpitopeFusion {
    sv_ids: [i32; 2],
    chromosomes: [String; 2],
    positions: [i32; 2],
}

pub struct NeoEpitopeWriter {
    output_dir: String,
    writer: Option<std::io::BufWriter<std::fs::File>>,
    gene_cache: &'static GeneCache,
    known_fusion_cache: &'static KnownFusionCache,
    fusions: Vec<NeoEpitopeFusion>,
}

const HEADER: &str = "geneIdUp\tgeneNameUp\tchromosomeUp\tpositionUp\torientationUp\tsvIdUp\tgeneIdDown\tgeneNameDown\tchromosomeDown\tpositionDown\torientationDown\tsvIdDown\tjunctionCopyNumber\tcopyNumber\tinsertSeq\tchainLength\ttranscriptsUp\ttranscriptsDown";

impl NeoEpitopeWriter {
    pub fn new(output_dir: &str, gene_cache: &'static GeneCache, known_fusion_cache: &'static KnownFusionCache) -> Self {
        NeoEpitopeWriter { output_dir: output_dir.to_string(), writer: None, gene_cache, known_fusion_cache, fusions: Vec::new() }
    }

    pub fn initialise_sample(&mut self, sample_id: &str) -> Result<()> {
        self.close()?;
        self.fusions.clear();
        let path = format!("{}/{}.linx.neoepitope.tsv", self.output_dir.trim_end_matches('/'), sample_id);
        let f = std::fs::File::create(&path).with_context(|| format!("error initialising neo-epitope output file: {path}"))?;
        let mut w = std::io::BufWriter::new(f);
        writeln!(w, "{HEADER}")?;
        self.writer = Some(w);
        Ok(())
    }

    pub fn close(&mut self) -> Result<()> {
        if let Some(mut w) = self.writer.take() {
            w.flush()?;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn process_fusion_candidate(
        &mut self, genes1: &[BreakendGene], genes2: &[BreakendGene], traversed: Option<&[Pair]>, disruption_finder: Option<&DisruptionFinder>,
        lower_link: Option<&Pair>, upper_link: Option<&Pair>, sv_copy_number: f64,
    ) -> Result<()> {
        if genes1.is_empty() || genes2.is_empty() {
            return Ok(());
        }
        let mut extension = [0i32; 2];
        let mut chain_length = 0;
        if let Some(t) = traversed {
            chain_length = t.iter().map(|x| x.base_length()).sum::<i32>();
            if t.len() as i32 > FUSION_MAX_CHAIN_LINKS || chain_length as i64 > FUSION_MAX_CHAIN_LENGTH {
                return Ok(());
            }
        }
        let mut fusion_added = false;
        for g1 in genes1 {
            for g2 in genes2 {
                if g1.is_upstream() == g2.is_upstream() {
                    continue;
                }
                let up = if g1.is_upstream() { g1 } else { g2 };
                if !is_candidate_upstream_gene(up) {
                    continue;
                }
                let down = if up == g1 { g2 } else { g1 };
                if !has_valid_traversal(up, traversed, g1.is_upstream(), disruption_finder) {
                    continue;
                }
                let gi = if g1.is_upstream() { FS_UP } else { FS_DOWN };
                extension[gi] = lower_link.map(|l| l.position_distance()).unwrap_or(-1);
                extension[switch_stream(gi)] = upper_link.map(|l| l.position_distance()).unwrap_or(-1);
                let max_upstream = self.max_upstream_distance(up.gene_name(), down.gene_name());
                let valid_up = self.find_valid_transcripts(up, extension[FS_UP], 0);
                let valid_down = self.find_valid_transcripts(down, extension[FS_DOWN], max_upstream);
                if valid_up.is_empty() || valid_down.is_empty() {
                    continue;
                }
                if self.is_duplicate(&genes1[0], &genes2[0]) || self.is_duplicate(&genes2[0], &genes1[0]) {
                    continue;
                }
                let avg_jcn = (up.jcn() + down.jcn()) * 0.5;
                let trans_up: Vec<&str> = valid_up.iter().map(|x| x.trans_name()).collect();
                let trans_down: Vec<&str> = valid_down.iter().map(|x| x.trans_name()).collect();
                let line = [
                    up.gene_id().to_string(),
                    up.gene_name().to_string(),
                    up.chromosome(),
                    up.position().to_string(),
                    up.orientation().to_string(),
                    up.var_id().to_string(),
                    down.gene_id().to_string(),
                    down.gene_name().to_string(),
                    down.chromosome(),
                    down.position().to_string(),
                    down.orientation().to_string(),
                    down.var_id().to_string(),
                    java_fixed(avg_jcn, 4),
                    java_fixed(sv_copy_number, 4),
                    up.insert_sequence(),
                    chain_length.to_string(),
                    trans_up.join(";"),
                    trans_down.join(";"),
                ]
                .join("\t");
                if let Some(w) = self.writer.as_mut() {
                    writeln!(w, "{line}")?;
                }
                if !fusion_added {
                    self.fusions.push(NeoEpitopeFusion {
                        sv_ids: [up.var_id(), down.var_id()],
                        chromosomes: [up.chromosome(), down.chromosome()],
                        positions: [up.position(), down.position()],
                    });
                    fusion_added = true;
                }
            }
        }
        Ok(())
    }

    fn find_valid_transcripts(&self, gene: &BreakendGene, link_extension: i32, max_upstream: i32) -> Vec<BreakendTrans> {
        let mut out = Vec::new();
        for t in gene.transcripts() {
            if gene.is_upstream() {
                if !t.is_disruptive() || t.trans_data.coding_start.is_none() {
                    continue;
                }
            } else {
                if t.trans_data.exons.len() <= 1 {
                    continue;
                }
                if t.region_type() == TranscriptRegionType::EXONIC && t.exon_downstream == t.exon_count() {
                    continue;
                }
                if t.trans_data.bio_type == BIOTYPE_NONSENSE_MED_DECAY {
                    continue;
                }
                if t.distance_upstream() > max_upstream {
                    continue;
                }
                let pre_trans_distance =
                    if gene.strand() == POS_STRAND { t.trans_start() - gene.position() } else { gene.position() - t.trans_end() };
                if pre_trans_distance > 0 {
                    let sa_pos = self.gene_cache.ensembl.find_preceding_gene_splice_acceptor_position(t.trans_id());
                    if sa_pos > 0 {
                        let sa_distance = if gene.strand() == POS_STRAND { t.trans_start() - sa_pos } else { sa_pos - t.trans_end() };
                        if sa_distance < pre_trans_distance {
                            continue;
                        }
                    }
                }
            }
            if link_extension > 0 {
                let trans_length = if (gene.strand() == POS_STRAND) == gene.is_upstream() {
                    (gene.position() - t.trans_data.trans_start).abs()
                } else {
                    (gene.position() - t.trans_data.trans_end).abs()
                };
                if link_extension < trans_length {
                    continue;
                }
            }
            out.push(t);
        }
        out
    }

    fn max_upstream_distance(&self, up_gene: &str, down_gene: &str) -> i32 {
        let k = self.known_fusion_cache;
        if k.has_known_fusion(up_gene, down_gene) {
            return MAX_UPSTREAM_DISTANCE_KNOWN;
        }
        if k.has_promiscuous_five_gene(up_gene) && k.is_high_impact_promiscuous(PROMISCUOUS_5, up_gene, down_gene) {
            return MAX_UPSTREAM_DISTANCE_KNOWN;
        }
        if k.has_promiscuous_three_gene(down_gene) && k.is_high_impact_promiscuous(PROMISCUOUS_3, up_gene, down_gene) {
            return MAX_UPSTREAM_DISTANCE_KNOWN;
        }
        MAX_UPSTREAM_DISTANCE_OTHER
    }

    fn is_duplicate(&self, g1: &BreakendGene, g2: &BreakendGene) -> bool {
        let c1 = g1.chromosome();
        let c2 = g2.chromosome();
        self.fusions.iter().any(|f| {
            [FS_UP, FS_DOWN].iter().any(|&fs| {
                let o = switch_stream(fs);
                f.sv_ids[fs] == g1.var_id()
                    && f.sv_ids[o] == g2.var_id()
                    && f.chromosomes[fs] == c1
                    && f.chromosomes[o] == c2
                    && f.positions[fs] == g1.position()
                    && f.positions[o] == g2.position()
            })
        })
    }
}

fn has_valid_traversal(up: &BreakendGene, traversed: Option<&[Pair]>, lower_to_upper: bool, finder: Option<&DisruptionFinder>) -> bool {
    let Some(traversed) = traversed else { return true };
    if traversed.is_empty() {
        return true;
    }
    let strand = up.strand() as i32;
    let finder = finder.expect("disruption finder for chained candidates");
    for pair in traversed {
        let orient = if lower_to_upper { pair.first_breakend().orientation() } else { pair.second_breakend().orientation() } as i32;
        let direction = if orient != strand { strand } else { -strand };
        if finder.pair_traverses_gene(pair, direction, false) {
            return false;
        }
    }
    true
}

fn is_candidate_upstream_gene(gene: &BreakendGene) -> bool {
    gene.transcripts().iter().any(|x| x.is_coding() && x.is_disruptive())
}
