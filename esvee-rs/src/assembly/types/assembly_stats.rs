//! AssemblyStats: read type counts and quality totals of an assembly.

use std::collections::HashMap;

use crate::assembly::constants::*;
use crate::assembly::read::read::Read;
use crate::assembly::read::read_utils::avg_base_quality_read;
use crate::assembly::types::junction::Junction;
use crate::assembly::types::support_read::SupportRead;
use crate::assembly::types::SupportType;
use crate::common::fmt::java_fixed;

/// Math.round(double)
pub fn java_round(v: f64) -> i64 {
    if v.is_nan() {
        return 0;
    }
    (v + 0.5).floor() as i64
}

#[derive(Default)]
pub struct AssemblyStats {
    pub junc_supps: i32,
    pub junc_mate_concordant: i32,
    pub junc_mate_discordant_remote: i32,
    pub junc_mate_discordant_ref_side: i32,
    pub junc_mate_ref_side: i32,
    pub junc_mate_unmapped_remote: i32,
    pub junc_mate_unmapped_ref_side: i32,
    pub indel_reads: i32,
    pub proximate_junc_read_ratio: f64,
    pub indel_length_total: i32,
    pub base_qual_total: i32,
    pub map_qual_total: i32,
    pub base_trim_total: i32,
    pub soft_clip_match_total: i32,
    pub soft_clip_mismatch_total: i32,
    pub ref_base_mismatch_total: i32,
    pub soft_clip_second_max_length: i32,
    pub max_ext_base_match_count: i32,
    pub non_extension_support_reads: i32,
    pub candidate_support_count: i32,
    pub unmapped_read_count: i32,
    pub read_count: i32,
    pub remote_supp_region_frequency: HashMap<String, i32>,
}

impl AssemblyStats {
    pub fn new() -> Self {
        AssemblyStats::default()
    }

    pub fn add_read(&mut self, support_read: &SupportRead, junction: &Junction, read: Option<&Read>) {
        self.read_count += 1;
        if support_read.type_().is_split_support() {
            if support_read.type_() == SupportType::Indel
                || (junction.indel_based() && support_read.indel_coords().is_some_and(|c| c.length >= MIN_VARIANT_LENGTH))
            {
                self.indel_reads += 1;
            }
            if support_read.is_supplementary() {
                self.junc_supps += 1;
            } else if support_read.is_paired_read() {
                let mate_past_junction = support_read.orientation().is_forward() == junction.is_forward();
                if support_read.is_mate_unmapped() {
                    if mate_past_junction {
                        self.junc_mate_unmapped_remote += 1;
                    } else {
                        self.junc_mate_unmapped_ref_side += 1;
                    }
                } else if support_read.is_discordant() {
                    if mate_past_junction {
                        self.junc_mate_discordant_remote += 1;
                    } else {
                        self.junc_mate_discordant_ref_side += 1;
                    }
                } else if (junction.is_forward() && support_read.mate_alignment_start() > junction.position)
                    || (!junction.is_forward() && support_read.mate_alignment_end() < junction.position)
                {
                    self.junc_mate_concordant += 1;
                } else {
                    self.junc_mate_ref_side += 1;
                }
            }
            self.soft_clip_match_total += support_read.extension_base_matches();
            self.soft_clip_mismatch_total += support_read.extension_base_mismatches();
            self.ref_base_mismatch_total += support_read.reference_mismatches().max(0);
            if support_read.extension_base_matches() > self.max_ext_base_match_count {
                self.soft_clip_second_max_length = self.max_ext_base_match_count;
                self.max_ext_base_match_count = support_read.extension_base_matches();
            } else if support_read.extension_base_matches() > self.soft_clip_second_max_length {
                self.soft_clip_second_max_length = support_read.extension_base_matches();
            }
            if junction.soft_clip_based() {
                self.track_supp_remote_regions(support_read, junction);
            }
        }
        self.map_qual_total += support_read.map_qual();
        self.base_trim_total += support_read.trim_count();
        if let Some(read) = read {
            self.base_qual_total += avg_base_quality_read(read);
            self.indel_length_total += read.cigar_elements().iter().filter(|x| x.op.is_indel()).map(|x| x.len).sum::<i32>();
        }
    }

    fn track_supp_remote_regions(&mut self, read: &SupportRead, junction: &Junction) {
        let Some(supp) = read.supplementary_data() else { return };
        if supp.map_quality < 30 {
            return;
        }
        if read.junction_extension_length(junction.orient) < ASSEMBLY_MIN_SOFT_CLIP_LENGTH {
            return;
        }
        if junction.is_forward() && read.alignment_end() != junction.position {
            return;
        } else if junction.is_reverse() && read.alignment_start() != junction.position {
            return;
        }
        let key = format!("{}_{}", supp.chromosome, java_round(supp.position as f64 / 1000.0));
        *self.remote_supp_region_frequency.entry(key).or_insert(0) += 1;
    }

    pub fn supp_remote_region_ratio(&self) -> f64 {
        if self.remote_supp_region_frequency.is_empty() {
            return 0.0;
        }
        if self.remote_supp_region_frequency.len() == 1 {
            return 1.0;
        }
        let mut total = 0;
        let mut max_region = 0;
        for count in self.remote_supp_region_frequency.values() {
            total += count;
            max_region = max_region.max(*count);
        }
        if total > 0 { max_region as f64 / total as f64 } else { 0.0 }
    }

    pub fn add_read_type_header(sj: &mut Vec<String>) {
        for h in ["JuncSupps", "JuncIndels", "JuncMateConcord", "JuncMateRefSide", "JuncMateDiscordRemote", "JuncMateDiscordRefSide", "JuncMateUnmappedRemote", "JuncMateUnmappedRefSide"] {
            sj.push(h.to_string());
        }
    }
    pub fn add_read_type_counts(&self, sj: &mut Vec<String>) {
        for v in [
            self.junc_supps,
            self.indel_reads,
            self.junc_mate_concordant,
            self.junc_mate_ref_side,
            self.junc_mate_discordant_remote,
            self.junc_mate_discordant_ref_side,
            self.junc_mate_unmapped_remote,
            self.junc_mate_unmapped_ref_side,
        ] {
            sj.push(v.to_string());
        }
    }
    pub fn add_read_stats_header(sj: &mut Vec<String>) {
        for h in [
            "SoftClipMatches", "SoftClipMismatches", "SoftClip2ndMaxLength", "RefBaseMismatches", "BaseTrimCount", "NonExtSupportReads",
            "ProxJuncReadRatio", "SuppRemoteLocationRatio", "AvgIndelLength", "AvgBaseQual", "AvgMapQual",
        ] {
            sj.push(h.to_string());
        }
    }
    pub fn add_read_stats(&self, sj: &mut Vec<String>) {
        sj.push(self.soft_clip_match_total.to_string());
        sj.push(self.soft_clip_mismatch_total.to_string());
        sj.push(self.soft_clip_second_max_length.to_string());
        sj.push(self.ref_base_mismatch_total.to_string());
        sj.push(self.base_trim_total.to_string());
        sj.push(self.non_extension_support_reads.to_string());
        sj.push(java_fixed(self.proximate_junc_read_ratio, 2));
        sj.push(java_fixed(self.supp_remote_region_ratio(), 2));
        sj.push(stat_string(self.indel_length_total, self.read_count as f64));
        sj.push(stat_string(self.base_qual_total, self.read_count as f64));
        sj.push(stat_string(self.map_qual_total, self.read_count as f64));
    }
}

fn stat_string(count: i32, read_count: f64) -> String {
    java_round(count as f64 / read_count).to_string()
}
