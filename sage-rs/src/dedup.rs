//! Variant de-duplication: DedupSnvMnv, DedupMixedGermlineSomatic, IndelDeduper, DedupMatching.

use crate::constants::*;
use crate::filters::{SageVariant, SoftFilter};
use crate::refdata::TranscriptData;
use crate::refseq::RefGenome;
use crate::types::*;
use std::collections::HashSet;

fn is_germline_and_not_tumor_filtered(filters: &HashSet<SoftFilter>) -> bool {
    !filters.is_empty() && filters.iter().all(|f| f.is_germline())
}

/// VariantDeduper.longerContainsShorter
fn longer_contains_shorter(shorter: &SimpleVariant, longer: &SimpleVariant) -> bool {
    if shorter.pos < longer.pos || shorter.end() > longer.end() {
        return false;
    }
    let offset = (shorter.pos - longer.pos) as usize;
    longer.alt.get(offset..offset + shorter.alt.len()).map_or(false, |s| s == &shorter.alt[..])
}

fn mnv_changed_bases(v: &SimpleVariant) -> usize {
    v.ref_.iter().zip(v.alt.iter()).filter(|(a, b)| a != b).count()
}

/// dedupMnvOverlaps
fn dedup_mnv_overlaps(variants: &mut [SageVariant]) {
    let mnvs: Vec<usize> = (0..variants.len())
        .filter(|&i| variants[i].variant().is_mnv() && variants[i].is_passing() && variants[i].has_local_phase_sets())
        .collect();
    let mut index = 0;
    while index + 1 < mnvs.len() {
        let first = mnvs[index];
        let max_pos = variants[first].position() + variants[first].variant().alt_len() - 1;
        let mut next_index = index + 1;
        while next_index < mnvs.len() {
            let next = mnvs[next_index];
            if variants[next].position() > max_pos {
                break;
            }
            if !variants[next].has_matching_lps(variants[first].local_phase_sets()) {
                next_index += 1;
                continue;
            }
            // dedupMnvPair
            let (fc, sc) = (mnv_changed_bases(variants[first].variant()), mnv_changed_bases(variants[next].variant()));
            let (fl, sl) = (variants[first].variant().ref_len(), variants[next].variant().ref_len());
            let filter_second = if fc != sc {
                fc > sc
            } else if fl != sl {
                fl < sl
            } else {
                variants[first].total_quality() > variants[next].total_quality()
            };
            let target = if filter_second { next } else { first };
            variants[target].filters.insert(SoftFilter::DedupMnv);
            if !variants[first].is_passing() {
                break;
            }
            next_index += 1;
        }
        index += 1;
    }
}

/// dedupMnvSnvs
fn dedup_mnv_snvs(variants: &mut [SageVariant]) {
    let candidates: Vec<usize> = (0..variants.len())
        .filter(|&i| !variants[i].is_indel() && variants[i].is_passing() && variants[i].has_local_phase_sets())
        .collect();
    let mut index = 0;
    while index + 1 < candidates.len() {
        let first = candidates[index];
        let max_pos = variants[first].position() + variants[first].variant().alt_len() - 1;
        let mut next_index = index + 1;
        while next_index < candidates.len() {
            let next = candidates[next_index];
            if variants[next].position() > max_pos {
                break;
            }
            if !variants[next].has_matching_lps(variants[first].local_phase_sets()) {
                next_index += 1;
                continue;
            }
            let first_snv = variants[first].variant().is_snv();
            if first_snv != variants[next].variant().is_snv() {
                let (snv, mnv) = if first_snv { (first, next) } else { (next, first) };
                if longer_contains_shorter(variants[snv].variant(), variants[mnv].variant()) {
                    variants[snv].filters.insert(SoftFilter::DedupSnvMnv);
                }
            }
            if !variants[first].is_passing() {
                break;
            }
            next_index += 1;
        }
        index += 1;
    }
}

/// TranscriptUtils.calcExonicCodingPhase
fn calc_exonic_coding_phase(exon_start: i32, exon_end: i32, phase_start: i32, coding_start: i32, coding_end: i32, strand: i8, position: i32) -> i32 {
    if !position_within(position, coding_start, coding_end) {
        return -1;
    }
    let pos_strand = strand == 1;
    let coding_bases =
        if pos_strand { position - exon_start.max(coding_start) + 1 } else { exon_end.min(coding_end) - position + 1 };
    if (pos_strand && coding_start == exon_start) || (!pos_strand && coding_end == exon_end) {
        let start_phase = if phase_start == -1 { 0 } else { phase_start };
        (start_phase + coding_bases) % 3
    } else if (pos_strand && coding_start < exon_start) || (!pos_strand && coding_end > exon_end) {
        (phase_start + coding_bases) % 3
    } else {
        coding_bases % 3
    }
}

fn codon_differences(codon: BaseRegion, v: &SimpleVariant) -> i32 {
    if codon.start > v.end() || v.pos > codon.end {
        return 0;
    }
    let s = codon.start.max(v.pos);
    let e = codon.end.min(v.end());
    (s..=e).filter(|p| v.ref_[(p - v.pos) as usize] != v.alt[(p - v.pos) as usize]).count() as i32
}

struct MixedGermlineDedup<'a> {
    transcripts: &'a [TranscriptData],
    index: usize,
    next_mixed_id: i32,
}

impl<'a> MixedGermlineDedup<'a> {
    /// TranscriptSelector.select + findCodon
    fn find_codon(&mut self, position: i32) -> Option<BaseRegion> {
        if self.transcripts.is_empty() {
            return None;
        }
        while self.index > 0 && self.transcripts[self.index].trans_start > position {
            self.index -= 1;
        }
        while self.index < self.transcripts.len() - 1 && self.transcripts[self.index].trans_end < position {
            self.index += 1;
        }
        let t = &self.transcripts[self.index];
        if !position_within(position, t.trans_start, t.trans_end) {
            return None;
        }
        let (cs, ce) = (t.coding_start?, t.coding_end?);
        if !position_within(position, cs, ce) {
            return None;
        }
        let exon = t.exons.iter().find(|x| position_within(position, x.start, x.end))?;
        let phase = calc_exonic_coding_phase(exon.start, exon.end, exon.phase_start, cs, ce, t.strand, position);
        let (mut start, mut end) = (position, position);
        if t.strand == 1 {
            match phase {
                1 => end = exon.end.min(position + 2),
                2 => {
                    start = exon.start.max(position - 1);
                    end = exon.end.min(position + 1);
                }
                0 => start = exon.start.max(position - 2),
                _ => {}
            }
        } else {
            match phase {
                1 => start = exon.start.max(position - 2),
                2 => {
                    start = exon.start.max(position - 1);
                    end = exon.end.min(position + 1);
                }
                0 => end = exon.end.min(position + 2),
                _ => {}
            }
        }
        Some(BaseRegion::new(start, end))
    }

    fn dedup(&mut self, variants: &mut [SageVariant]) {
        let germline_snv = |v: &SageVariant| v.variant().is_snv() && !v.is_passing() && is_germline_and_not_tumor_filtered(&v.filters);
        let passing_mnv = |v: &SageVariant| v.variant().is_mnv() && v.is_passing() && v.has_local_phase_sets();
        let passing_snv = |v: &SageVariant| v.variant().is_snv() && v.is_passing() && v.has_local_phase_sets();
        let candidates: Vec<usize> =
            (0..variants.len()).filter(|&i| germline_snv(&variants[i]) || passing_snv(&variants[i]) || passing_mnv(&variants[i])).collect();
        let mut index = 0;
        while index + 1 < candidates.len() {
            let first = candidates[index];
            let max_pos = variants[first].position() + variants[first].variant().alt_len() - 1;
            let (mut mnv, mut snv_passing, mut snv_germline) = (None, None, None);
            if passing_mnv(&variants[first]) {
                mnv = Some(first);
            } else if passing_snv(&variants[first]) {
                snv_passing = Some(first);
            } else {
                snv_germline = Some(first);
            }
            let mut next_index = index + 1;
            while next_index < candidates.len() {
                let next = candidates[next_index];
                if variants[next].position() - max_pos > 5 {
                    break;
                }
                if mnv.is_none() && passing_mnv(&variants[next]) {
                    mnv = Some(next);
                } else if snv_passing.is_none() && passing_snv(&variants[next]) {
                    snv_passing = Some(next);
                } else if snv_germline.is_none() && germline_snv(&variants[next]) {
                    snv_germline = Some(next);
                }
                if let (Some(m), Some(sp), Some(sg)) = (mnv, snv_passing, snv_germline) {
                    self.check_group(variants, m, sp, sg);
                    break;
                }
                next_index += 1;
            }
            index += 1;
        }
    }

    fn check_group(&mut self, variants: &mut [SageVariant], mnv: usize, snv_passing: usize, snv_germline: usize) {
        if !variants[snv_passing].has_matching_lps(variants[mnv].local_phase_sets()) {
            return;
        }
        if !longer_contains_shorter(variants[snv_germline].variant(), variants[mnv].variant()) {
            return;
        }
        if variants[mnv].mixed_impact == 0 {
            variants[mnv].mixed_impact = self.next_mixed_id;
            self.next_mixed_id += 1;
        }
        variants[snv_germline].mixed_impact = variants[mnv].mixed_impact;
        if longer_contains_shorter(variants[snv_passing].variant(), variants[mnv].variant()) {
            let codon = self.find_codon(variants[snv_passing].position());
            let keep_mnv = codon.map_or(false, |c| {
                codon_differences(c, variants[mnv].variant()) > codon_differences(c, variants[snv_passing].variant())
            });
            let target = if keep_mnv { snv_passing } else { mnv };
            variants[target].filters.insert(SoftFilter::DedupMixedGermlineSomatic);
        }
    }
}

struct VariantData {
    /// index into the variant list
    v: usize,
    indel_score: i32,
    core_pos_end: i32,
    flank_pos_start: i32,
    flank_pos_end: i32,
}

fn calc_indel_score(v: &SageVariant) -> i32 {
    let rc = &v.tumor_counters[0];
    let edge = jround(rc.read_edge_distance.max_alt_distance_from_edge * DEFAULT_READ_LENGTH as f64) as i32;
    (if v.is_indel() { v.variant().indel_length_abs() * 8 } else { 0 }) + edge - 12 * rc.min_number_of_events
}

fn allow_by_filter(v: &SageVariant) -> bool {
    v.is_passing() || is_germline_and_not_tumor_filtered(&v.filters)
}

fn read_edge_distance_threshold(tier: Tier) -> f64 {
    if matches!(tier, Tier::Hotspot | Tier::Panel) {
        0.2
    } else {
        0.33
    }
}

fn is_dedup_candidate(variants: &[SageVariant], indel: &VariantData, d: &VariantData, require_core_end: bool) -> bool {
    let v = &variants[d.v];
    if require_core_end {
        if positions_within(v.position(), d.core_pos_end, indel.flank_pos_start, indel.flank_pos_end) {
            return true;
        }
    } else if position_within(v.position(), indel.flank_pos_start, indel.flank_pos_end) {
        return true;
    }
    v.tumor_counters[0].read_edge_distance.max_alt_distance_from_edge < read_edge_distance_threshold(v.tier())
}

/// IndelDeduper.buildAltBasesString
pub fn build_alt_bases(ref_bases: &[u8], ref_pos_start: i32, ref_pos_end: i32, variants: &[&SimpleVariant]) -> Vec<u8> {
    let mut sorted: Vec<&SimpleVariant> = variants.to_vec();
    sorted.sort_by(|a, b| b.pos.cmp(&a.pos)); // stable, descending position
    let mut alt = ref_bases.to_vec();
    for v in sorted {
        if !position_within(v.pos, ref_pos_start, ref_pos_end) {
            continue;
        }
        let rel = ((v.pos - ref_pos_start) as usize).min(alt.len());
        let mut new_alt = alt[..rel].to_vec();
        new_alt.extend_from_slice(&v.alt);
        let gap = if v.is_insert() { 1 } else { v.ref_.len() };
        if rel + gap <= alt.len() {
            new_alt.extend_from_slice(&alt[rel + gap..]);
        }
        alt = new_alt;
    }
    alt
}

fn contains(hay: &[u8], needle: &[u8]) -> bool {
    needle.is_empty() || hay.windows(needle.len()).any(|w| w == needle)
}

struct IndelDeduper<'a> {
    genome: &'a RefGenome,
    chr: &'a str,
    iterations: i32,
}

impl<'a> IndelDeduper<'a> {
    fn dedup(&mut self, variants: &mut [SageVariant]) {
        let mut data: Vec<VariantData> = Vec::new();
        let mut indels: Vec<usize> = Vec::new();
        let mut candidates: Vec<usize> = Vec::new();
        for (i, v) in variants.iter().enumerate() {
            if !v.has_local_phase_sets() {
                continue;
            }
            let ctx = &v.tumor_counters[0].ctx;
            data.push(VariantData {
                v: i,
                indel_score: if v.is_indel() { calc_indel_score(v) } else { 0 },
                core_pos_end: ctx.core_position_end,
                flank_pos_start: ctx.alignment_start,
                flank_pos_end: ctx.alignment_end,
            });
            let d = data.len() - 1;
            if v.is_indel() && allow_by_filter(v) {
                indels.push(d);
            }
            candidates.push(d);
        }
        if indels.is_empty() {
            return;
        }
        indels.sort_by(|&a, &b| data[b].indel_score.cmp(&data[a].indel_score));
        candidates.sort_by(|&a, &b| data[a].flank_pos_start.cmp(&data[b].flank_pos_start));

        while !indels.is_empty() {
            let indel = indels.remove(0);
            // findDedupGroup
            let mut group: Vec<usize> = Vec::new();
            for &c in &candidates {
                if c == indel {
                    continue;
                }
                let v = &variants[data[c].v];
                if v.variant().position_end() < data[indel].flank_pos_start - 60 {
                    continue;
                }
                if v.position() > data[indel].flank_pos_end + 60 {
                    break;
                }
                if !variants[data[indel].v].has_matching_lps(v.local_phase_sets()) {
                    continue;
                }
                if is_dedup_candidate(variants, &data[indel], &data[c], false) {
                    group.push(c);
                }
            }
            if group.is_empty() {
                continue;
            }
            self.dedup_indel_group(variants, &data, indel, &mut group);
            for d in group {
                if let Some(p) = candidates.iter().position(|&x| x == d) {
                    candidates.remove(p);
                }
                if variants[data[d].v].is_indel() {
                    if let Some(p) = indels.iter().position(|&x| x == d) {
                        indels.remove(p);
                    }
                }
            }
        }
    }

    fn dedup_indel_group(&mut self, variants: &mut [SageVariant], data: &[VariantData], indel: usize, group: &mut Vec<usize>) {
        let ind = &data[indel];
        if group.iter().any(|&d| !variants[data[d].v].is_indel() && calc_indel_score(&variants[data[d].v]) > ind.indel_score) {
            return;
        }
        let indel_position = variants[ind.v].position();
        let mut overlapped: Vec<usize> = Vec::new();
        let mut has_valid = false;
        let mut index = 0;
        while index < group.len() {
            let v = &variants[data[group[index]].v];
            has_valid |= allow_by_filter(v);
            if v.variant().is_delete()
                && positions_overlap(indel_position, indel_position, v.position(), v.variant().position_end())
            {
                overlapped.push(group.remove(index));
            } else {
                index += 1;
            }
        }
        if !has_valid {
            return;
        }
        let Ok(ref_bases) = self.genome.bases(self.chr, ind.flank_pos_start, ind.flank_pos_end) else { return };
        if ref_bases.is_empty() {
            return;
        }
        let indel_bases = variants[ind.v].tumor_counters[0].ctx.read_bases.clone();
        self.iterations = 0;
        let mut deduped: Vec<usize> = Vec::new();
        if !self.check_combinations(variants, data, indel, group, &mut deduped, &indel_bases, &ref_bases) {
            for &d in &overlapped {
                variants[data[d].v].filters.insert(SoftFilter::DedupIndel);
            }
            group.extend(&overlapped);
            return;
        }
        group.extend(&overlapped);
        deduped.extend(&overlapped);
        group.push(indel);
        let non_deduped: Vec<usize> = group.iter().copied().filter(|d| !deduped.contains(d)).collect();
        for &d in group.iter() {
            let vi = data[d].v;
            let indel_passing = variants[ind.v].is_passing();
            let indel_filters: Vec<SoftFilter> = variants[ind.v].filters.iter().copied().collect();
            if deduped.contains(&d) {
                if !is_dedup_candidate(variants, ind, &data[d], true) {
                    continue;
                }
                variants[vi].filters.insert(SoftFilter::DedupIndel);
                if !indel_passing {
                    variants[vi].filters.extend(indel_filters);
                }
            } else if is_dedup_candidate(variants, ind, &data[d], true) && !indel_passing {
                if let Some(r) = variants[ind.v].ref_counters.first() {
                    let ref_full = r.counts.full;
                    let ref_mismatch = r.simple_alt_matches + r.counts.core;
                    if ref_full > ref_mismatch {
                        variants[vi].filters.extend(indel_filters);
                    }
                }
            } else if indel_passing && recover_filtered_variant(variants, data, vi, &non_deduped) {
                variants[vi].filters.clear();
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn check_combination(
        &mut self, variants: &[SageVariant], data: &[VariantData], indel: usize, all: &[usize], selected: &[usize],
        deduped: &mut Vec<usize>, indel_bases: &[u8], ref_bases: &[u8],
    ) -> bool {
        self.iterations += 1;
        let ind = &data[indel];
        let mut test: Vec<&SimpleVariant> = vec![variants[ind.v].variant()];
        test.extend(selected.iter().map(|&d| variants[data[d].v].variant()));
        let net = build_alt_bases(ref_bases, ind.flank_pos_start, ind.flank_pos_end, &test);
        let matched = if indel_bases.len() >= net.len() { contains(indel_bases, &net) } else { contains(&net, indel_bases) };
        if matched {
            deduped.extend(all.iter().filter(|d| !selected.contains(d)));
        }
        matched
    }

    #[allow(clippy::too_many_arguments)]
    fn check_combinations(
        &mut self, variants: &[SageVariant], data: &[VariantData], indel: usize, group: &[usize], deduped: &mut Vec<usize>,
        indel_bases: &[u8], ref_bases: &[u8],
    ) -> bool {
        if self.check_combination(variants, data, indel, group, &[], deduped, indel_bases, ref_bases) {
            return true;
        }
        // VariantData.compareTo
        let cmp = |a: &usize, b: &usize| -> std::cmp::Ordering {
            let (va, vb) = (&variants[data[*a].v], &variants[data[*b].v]);
            va.position()
                .cmp(&vb.position())
                .then(va.variant().ref_len().cmp(&vb.variant().ref_len()))
                .then(vb.tumor_counters[0].alt_support().cmp(&va.tumor_counters[0].alt_support()))
        };
        let n = group.len();
        let max_selection = if n <= 6 { n } else { 3 };
        for size in 1..=max_selection {
            // Java enumerates subsets of a HashSet keyed by identity hash, so ties are ordered
            // arbitrarily there; here subsets follow the group order.
            let mut lists: Vec<Vec<usize>> = Vec::new();
            let mut mask: u64 = (1u64 << size) - 1;
            while mask < (1u64 << n) {
                let mut subset: Vec<usize> = (0..n).filter(|i| mask & (1 << i) != 0).map(|i| group[i]).collect();
                subset.sort_by(cmp);
                lists.push(subset);
                // next larger integer with the same number of set bits
                let c = mask & mask.wrapping_neg();
                let r = mask + c;
                mask = (((r ^ mask) >> 2) / c) | r;
            }
            lists.sort_by(|a, b| {
                for i in 0..a.len() {
                    let o = cmp(&a[i], &b[i]);
                    if o != std::cmp::Ordering::Equal {
                        return o;
                    }
                }
                std::cmp::Ordering::Equal
            });
            for list in &lists {
                if self.check_combination(variants, data, indel, group, list, deduped, indel_bases, ref_bases) {
                    return true;
                }
            }
            if self.iterations >= 100 {
                deduped.extend(group);
                return true;
            }
        }
        false
    }
}

fn recover_filtered_variant(variants: &[SageVariant], data: &[VariantData], vi: usize, non_deduped: &[usize]) -> bool {
    let v = &variants[vi];
    if v.is_passing() {
        return false;
    }
    if v.is_indel() && v.local_phase_sets().len() == 1 && non_deduped.len() > 1 {
        let mut max_percent = 0.0;
        let lps = v.local_phase_sets()[0];
        for &d in non_deduped {
            if data[d].v == vi {
                continue;
            }
            let o = &variants[data[d].v].tumor_counters[0];
            let mut total = 0;
            let mut matched = 0;
            for i in 0..o.local_phase_sets.len() {
                total += o.lps_counts[i];
                if o.local_phase_sets[i] == lps {
                    matched = o.lps_counts[i];
                }
            }
            max_percent = jmax(max_percent, matched as f64 / total as f64);
        }
        if max_percent < 0.1 {
            return false;
        }
    }
    !v.filters.iter().any(|f| f.is_germline())
}

/// DedupMatching.dedupMatchingVariants
fn dedup_matching(variants: &mut [SageVariant]) {
    let passing: Vec<usize> = (0..variants.len()).filter(|&i| variants[i].is_passing()).collect();
    let mut index = 0;
    while index + 1 < passing.len() {
        let (first, next) = (passing[index], passing[index + 1]);
        if variants[first].variant() == variants[next].variant() {
            let (source, dest) =
                if variants[first].total_quality() >= variants[next].total_quality() { (next, first) } else { (first, next) };
            let (sets, counts) = {
                let t = &variants[source].tumor_counters[0];
                (t.local_phase_sets.clone(), t.lps_counts.clone())
            };
            for (lps, count) in sets.into_iter().zip(counts) {
                variants[dest].tumor_counters[0].add_local_phase_set(lps, count, 0.0);
            }
            variants[source].filters.insert(SoftFilter::DedupMatch);
            index += 2;
        } else {
            index += 1;
        }
    }
}

/// VariantDeduper.processVariants
pub fn process_variants(variants: &mut [SageVariant], transcripts: &[TranscriptData], genome: &RefGenome, chr: &str) {
    dedup_mnv_overlaps(variants);
    MixedGermlineDedup { transcripts, index: 0, next_mixed_id: 1 }.dedup(variants);
    dedup_mnv_snvs(variants);
    IndelDeduper { genome, chr, iterations: 0 }.dedup(variants);
    dedup_matching(variants);
}
