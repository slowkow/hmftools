//! JunctionTracker and JunctionUtils: the read groups of a partition, the junctions their reads establish and
//! the fragments that support each junction.

use super::depth::DepthTracker;
use super::discordant::{DiscordantGroups, DiscordantStats};
use super::filters::*;
use super::hotspot::KnownHotspot;
use super::jmap::JavaHashMap;
use super::rec::*;
use super::saga::SagaLocationMatcher;
use super::types::*;
use crate::common::compat;
use std::collections::HashSet;
use std::sync::Arc;

pub const INVALID_JUNC_INDEX: i32 = -1;
pub const SIMPLE_SEARCH_COUNT: usize = 10;
pub const SBX_MAX_DUPLICATE_DISTANCE: i32 = 2;

/// ReadIdTrimmer: drops the part of a read id before its third ':' while every id has it in the same place.
///
/// Java (finding P12): each partition and the spanning-read cache have their own trimmer, whose prefix is
/// the part before the third ':' of the first name it is given; after that only the position of that ':'
/// is checked, and the trimmer switches itself off at a name whose ':' is elsewhere. With reads of several
/// sequencing runs in the one BAM, names of different flowcells whose prefixes have the same length are cut
/// to the same `lane:tile:x:y` id, and a switch-off part way leaves earlier mates under their cut ids and
/// later ones under their whole names. Corrected: one prefix for the whole run, that of the BAM's first
/// record (`trim_prefix`), and only names that start with it are cut.
pub struct ReadIdTrimmer {
    enabled: bool,
    trim_index: i32,
    /// corrected: the prefix to cut, the same in every partition
    prefix: Option<Arc<[u8]>>,
}

/// Corrected (finding P12): the part of a read name up to its third ':', the prefix for the whole run. None
/// for a name without three ':' or with nothing after the third, for which Java switches trimming off.
pub fn trim_prefix(read_name: &[u8]) -> Option<Arc<[u8]>> {
    let index = find_read_id_trim_index(read_name);
    (index >= 0 && (index as usize) < read_name.len()).then(|| Arc::from(&read_name[..index as usize]))
}

impl ReadIdTrimmer {
    /// Java's trimmer, which takes its prefix from the first name it is given.
    pub fn new(enabled: bool) -> Self {
        ReadIdTrimmer { enabled, trim_index: -1, prefix: None }
    }
    /// The run's trimmer: corrected, the prefix of the run (none: no trimming); with finding P12 reproduced,
    /// Java's.
    pub fn for_run(enabled: bool, prefix: Option<&Arc<[u8]>>) -> Self {
        if compat::java("P12") {
            return Self::new(enabled);
        }
        ReadIdTrimmer { enabled: enabled && prefix.is_some(), trim_index: -1, prefix: prefix.cloned() }
    }
    pub fn enabled(&self) -> bool {
        self.enabled
    }
    pub fn trim<'a>(&mut self, read_id: &'a [u8]) -> &'a [u8] {
        if !self.enabled {
            return read_id;
        }
        if let Some(prefix) = &self.prefix {
            return match read_id.strip_prefix(&prefix[..]) {
                Some(rest) if !rest.is_empty() => rest,
                _ => read_id,
            };
        }
        if self.trim_index < 0 {
            self.trim_index = find_read_id_trim_index(read_id);
            if self.trim_index < 0 {
                self.enabled = false;
                return read_id;
            }
        }
        let t = self.trim_index as usize;
        if t >= read_id.len() || read_id[t - 1] != b':' {
            self.enabled = false;
            return read_id;
        }
        &read_id[t..]
    }
}

fn find_read_id_trim_index(read_id: &[u8]) -> i32 {
    let mut delim_count = 0;
    for (i, &c) in read_id.iter().enumerate() {
        if c == b':' {
            delim_count += 1;
        }
        if delim_count == 3 {
            return i as i32 + 1;
        }
    }
    -1
}

/// The settings of a run that the tracker reads.
#[derive(Clone)]
pub struct TrackerConfig {
    pub filter: ReadFilterConfig,
    pub trim_read_id: bool,
    /// corrected (finding P12): the read name prefix of the run
    pub trim_prefix: Option<Arc<[u8]>>,
    pub contigs: Arc<Contigs>,
}

/// JunctionUtils.findJunctionIndex: the index of a junction at the position, else of the one before it, else
/// invalid. `position_at` gives the position of the junction at an index of the ordered list.
pub fn find_junction_index(count: usize, position_at: impl Fn(usize) -> i32, position: i32, simple_search_count: usize) -> i32 {
    if count == 0 || position < position_at(0) {
        return INVALID_JUNC_INDEX;
    }
    let max_index = count as i32 - 1;
    if position >= position_at(count - 1) {
        return max_index;
    }
    if count <= simple_search_count {
        for index in 0..count {
            let p = position_at(index);
            if p < position {
                continue;
            }
            if p > position {
                return index as i32 - 1;
            }
            return index as i32;
        }
        return max_index;
    }
    let mut current = (count / 2) as i32;
    let mut lower = 0i32;
    let mut upper = max_index;
    let mut iterations = 0;
    loop {
        let p = position_at(current as usize);
        if p == position {
            return current;
        }
        if upper == lower + 1 {
            break;
        }
        if position < p {
            if current == 0 {
                break;
            }
            if current == 1 {
                if lower == current {
                    break;
                }
                upper = 1;
                current = 0;
            } else {
                upper = current;
                current = (lower + upper) / 2;
            }
        } else {
            if current == max_index {
                break;
            }
            if current == max_index - 1 {
                if upper == current {
                    break;
                }
                lower = current;
                current = max_index;
            } else {
                lower = current;
                current = (lower + upper) / 2;
            }
        }
        iterations += 1;
        if iterations > 100 {
            crate::warn!("junction index search iterations({iterations}) junctions({count}) index(cur={current} low={lower} high={upper})");
            break;
        }
    }
    if position > position_at(current as usize) { current } else { current - 1 }
}

/// JunctionUtils.hasOtherJunctionSupport
pub fn has_other_junction_support(read: &PrepRead, junction: &JunctionData, filter: &ReadFilterConfig) -> bool {
    let (unclipped_start, unclipped_end) = (read.unclipped_start, read.unclipped_end);
    if position_within(junction.position, unclipped_start, unclipped_end) {
        let junction_distance = if junction.is_forward() {
            (unclipped_end - junction.position).abs().min((read.end - junction.position).abs())
        } else {
            (unclipped_start - junction.position).abs().min((read.start - junction.position).abs())
        };
        if junction.is_forward() && !read.is_right_clipped() {
            return false;
        } else if junction.is_reverse() && !read.is_left_clipped() {
            return false;
        }
        let max_permitted = if junction.internal_indel { MAX_JUNC_SC_VS_INDEL_POSITION_DIFF } else { MAX_JUNC_POSITION_DIFF };
        return junction_distance <= max_permitted;
    }
    if junction.orient != read.orientation() {
        return false;
    }
    let junction_distance;
    if junction.is_forward() {
        if read.end > junction.position {
            return false;
        }
        junction_distance = (read.end - junction.position).abs();
    } else {
        if read.start < junction.position {
            return false;
        }
        junction_distance = (read.start - junction.position).abs();
    }
    if junction_distance <= filter.max_supporting_fragment_distance {
        return is_chimeric_read(&read.rec, filter);
    }
    false
}

/// JunctionUtils.hasExactJunctionSupport: clipped at the junction base, or clipped close by with soft-clip and
/// reference bases that agree with those of the junction's top read
pub fn has_exact_junction_support(read: &PrepRead, junction: &JunctionData, junc_read: Option<&PrepRead>) -> bool {
    let (read_start, read_end) = (read.start, read.end);
    let read_length = read.bases.len() as i32;
    let Some(junc_read) = junc_read else { return false };
    let junc_read_length = junc_read.bases.len() as i32;
    if junc_read.cigar().len() < 2 {
        return false;
    }
    let (read_index_start, junc_read_index_start);
    if junction.is_forward() {
        if !read.is_right_clipped() {
            return false;
        }
        if read_end == junction.position {
            return true;
        }
        if (read_end - junction.position).abs() > MAX_JUNC_POSITION_DIFF {
            return false;
        }
        if !(read_start < junction.position && read.unclipped_end > junction.position) {
            return false;
        }
        let junc_pos_diff = junc_read.end - read.end;
        if junc_pos_diff < 0 {
            read_index_start = read_length - read.sc_right + junc_pos_diff;
            junc_read_index_start = junc_read_length - junc_read.sc_right;
        } else {
            read_index_start = read_length - read.sc_right;
            junc_read_index_start = junc_read_length - junc_read.sc_right - junc_pos_diff;
        }
    } else {
        if !read.is_left_clipped() {
            return false;
        }
        if read_start == junction.position {
            return true;
        }
        if (read_start - junction.position).abs() > MAX_JUNC_POSITION_DIFF {
            return false;
        }
        if !(read_end > junction.position && read.unclipped_start < junction.position) {
            return false;
        }
        let junc_pos_diff = junc_read.start - read.start;
        if junc_pos_diff < 0 {
            read_index_start = read.sc_left - 1;
            junc_read_index_start = junc_read.sc_left - 1 - junc_pos_diff;
        } else {
            read_index_start = read.sc_left - 1 + junc_pos_diff;
            junc_read_index_start = junc_read.sc_left - 1;
        }
    }
    bases_meet_exact_support(junc_read, read, junc_read_index_start, read_index_start, junction.is_reverse())
}

fn bases_meet_exact_support(junc_read: &PrepRead, read: &PrepRead, mut junc_read_index: i32, mut read_index: i32, search_down: bool) -> bool {
    let (mut high_qual_mismatches, mut base_matches, mut bases_compared) = (0, 0, 0);
    let read_index_max = read.bases.len() as i32 - 1;
    let junc_read_index_max = junc_read.bases.len() as i32 - 1;
    let (read_quals, junc_quals) = (read.base_qualities(), junc_read.base_qualities());
    let step = if search_down { -1 } else { 1 };
    let t = crate::common::base_qual::low_base_qual_threshold();
    while junc_read_index >= 0 && read_index >= 0 && junc_read_index <= junc_read_index_max && read_index <= read_index_max {
        bases_compared += 1;
        if read.bases[read_index as usize] == junc_read.bases[junc_read_index as usize] {
            base_matches += 1;
        } else if !below_qual_threshold(read_quals[read_index as usize], t) && !below_qual_threshold(junc_quals[junc_read_index as usize], t) {
            high_qual_mismatches += 1;
            if high_qual_mismatches > MAX_HIGH_QUAL_BASE_MISMATCHES {
                return false;
            }
        }
        junc_read_index += step;
        read_index += step;
    }
    // no bases compared gives a NaN, which fails the test as it does in Java
    base_matches as f64 / bases_compared as f64 >= MIN_EXACT_BASE_PERC
}

fn above_aligned_score_difference(read: &PrepRead, min_align_score: i32) -> bool {
    let Some(as_score) = aux_int(&read.rec, b"AS") else { return false };
    if as_score < min_align_score as i64 {
        return false;
    }
    aux_int(&read.rec, b"XS").is_some_and(|xs| as_score - xs >= MIN_ALIGNMENT_SCORE_DIFF as i64)
}

/// FragmentData: a read with a supplementary alignment, as compared when marking supplementary duplicates
struct FragmentData {
    read: ReadId,
    is_primary: bool,
    read_is_lower_vs_supp: bool,
    other_position: i32,
    other_chr: i32,
    other_orient: Orientation,
    is_consensus: bool,
    orient: Orientation,
    unclipped_start: i32,
    unclipped_end: i32,
}

impl FragmentData {
    fn new(id: ReadId, read: &PrepRead, contigs: &Contigs) -> FragmentData {
        let supp = read.supp.as_ref().expect("supplementary data");
        let read_chr = contigs.name(read.chr);
        let read_is_lower_vs_supp =
            if read_chr == supp.chr_name { read.start <= supp.position } else { chromosome_rank(read_chr) < chromosome_rank(&supp.chr_name) };
        let other_position = if read.is_paired() {
            match aux_str(&read.rec, b"MC") {
                Some(mc) => five_prime_unclipped_str(mate_start(&read.rec), mc, read.mate_orientation().is_forward()),
                None => read.mate_pos_start,
            }
        } else if read.orientation().is_forward() {
            read.unclipped_end
        } else {
            read.unclipped_start
        };
        FragmentData {
            read: id,
            is_primary: !read.is_supplementary(),
            read_is_lower_vs_supp,
            other_position,
            other_chr: if read.is_paired() { read.mate_chr } else { read.chr },
            other_orient: if read.is_paired() { read.mate_orientation() } else { read.orientation() },
            is_consensus: has_aux(&read.rec, b"CR"),
            orient: read.orientation(),
            unclipped_start: read.unclipped_start,
            unclipped_end: read.unclipped_end,
        }
    }
    fn matches(&self, o: &FragmentData) -> bool {
        self.other_chr == o.other_chr && self.other_orient == o.other_orient && self.other_position == o.other_position
    }
    fn within_position_range(&self, o: &FragmentData, permitted: i32) -> bool {
        if self.orient != o.orient {
            return false;
        }
        (self.unclipped_start - o.unclipped_start).abs() + (self.unclipped_end - o.unclipped_end).abs() <= permitted
    }
}

pub struct JunctionTracker {
    pub region: Region,
    pub cfg: TrackerConfig,
    known_hotspots: Vec<Arc<KnownHotspot>>,
    pub reads: Vec<PrepRead>,
    /// slots of `reads` whose read was dropped with its group, for reuse
    free_reads: Vec<ReadId>,
    pub groups: Vec<ReadGroup>,
    /// every junction made; `junctions` is the ordered list of those in use
    pub jstore: Vec<JunctionData>,
    group_map: JavaHashMap<Box<[u8]>, GroupId>,
    expected_read_ids: HashSet<Box<[u8]>>,
    pub expected_groups: Vec<GroupId>,
    /// groups with reads in another partition that support no junction here; a set in Java
    remote_candidate_groups: Vec<GroupId>,
    candidate_discordant_groups: Vec<GroupId>,
    pub junctions: Vec<JuncId>,
    /// discordant junctions at a remote region outside the partition, kept apart from `junctions` when
    /// finding P8 is corrected
    pub outside_junctions: Vec<JuncId>,
    recent_junction_index: i32,
    trimmer: ReadIdTrimmer,
    pub initial_supporting_frags: i32,
    pub discordant_stats: DiscordantStats,
    discordant_finder: DiscordantGroups,
    saga: Option<SagaLocationMatcher>,
    pub duplicate_groups: i32,
}

impl JunctionTracker {
    pub fn new(region: Region, cfg: TrackerConfig, known_hotspots: Vec<Arc<KnownHotspot>>, saga: Option<SagaLocationMatcher>) -> Self {
        let discordant_finder = DiscordantGroups::new(region.clone(), cfg.filter.observed_frag_length_max, known_hotspots.clone(), cfg.contigs.clone());
        let trimmer = ReadIdTrimmer::for_run(cfg.trim_read_id, cfg.trim_prefix.as_ref());
        JunctionTracker {
            region,
            cfg,
            known_hotspots,
            reads: Vec::new(),
            free_reads: Vec::new(),
            groups: Vec::new(),
            jstore: Vec::new(),
            group_map: JavaHashMap::new(),
            expected_read_ids: HashSet::new(),
            expected_groups: Vec::new(),
            remote_candidate_groups: Vec::new(),
            candidate_discordant_groups: Vec::new(),
            junctions: Vec::new(),
            outside_junctions: Vec::new(),
            recent_junction_index: INVALID_JUNC_INDEX,
            trimmer,
            initial_supporting_frags: 0,
            discordant_stats: DiscordantStats::default(),
            discordant_finder,
            saga,
            duplicate_groups: 0,
        }
    }

    #[cfg(test)]
    pub fn junction(&self, index: usize) -> &JunctionData {
        &self.jstore[self.junctions[index] as usize]
    }

    /// The group of a read id, as held by the map.
    #[cfg(test)]
    pub fn group_of(&self, read_id: &[u8]) -> Option<GroupId> {
        self.group_map.get(&read_id.to_vec().into_boxed_slice()).copied()
    }

    /// processRead: adds the read to the group of its id. Returns its index, or none for a read that is dropped.
    pub fn process_read(&mut self, read: PrepRead) -> Option<ReadId> {
        if !read.has_mate() && !read.has_supp_alignment() && read.read_type == ReadType::NoSupport {
            return None;
        }
        let read_id: Box<[u8]> = self.trimmer.trim(read.id()).into();
        let rid = match self.free_reads.pop() {
            Some(slot) => {
                self.reads[slot as usize] = read;
                slot
            }
            None => {
                self.reads.push(read);
                (self.reads.len() - 1) as ReadId
            }
        };
        match self.group_map.get(&read_id).copied() {
            None => {
                let gid = self.groups.len() as GroupId;
                self.groups.push(ReadGroup::new(read_id.clone(), rid));
                self.group_map.insert(read_id, gid);
            }
            Some(gid) => {
                self.groups[gid as usize].reads.push(rid);
                let g = &self.groups[gid as usize];
                // a complete pair without any support is of no further use
                if group_is_simple_complete(g, &self.reads) && group_all_no_support(g, &self.reads) {
                    self.group_map.remove(&read_id);
                    let ids = std::mem::take(&mut self.groups[gid as usize].reads);
                    for i in ids {
                        let _ = self.reads[i as usize].take();
                        self.free_reads.push(i);
                    }
                }
            }
        }
        Some(rid)
    }

    pub fn set_expected_reads(&mut self, expected: HashSet<Box<[u8]>>) {
        self.expected_read_ids.extend(expected);
    }

    pub fn assign_fragments(&mut self, depth: &DepthTracker) {
        self.assign_junction_fragments_and_support();
        self.filter_junctions(depth);
        self.find_discordant_groups();
        if self.saga.is_some() {
            self.match_junctions_to_saga();
        }
    }

    /// The groups of the map in Java's iteration order.
    fn map_order(&self) -> Vec<GroupId> {
        self.group_map.values().copied().collect()
    }

    /// The step of assignJunctionFragmentsAndSupport that takes the read ids another partition expects here:
    /// in map order, every group that is not a duplicate and is expected is marked as such. It reads and
    /// changes nothing the rest of the assignment uses, so a run applies it once the earlier partitions are
    /// complete.
    pub fn mark_expected_groups(&mut self) {
        if self.expected_read_ids.is_empty() {
            return;
        }
        for gid in self.map_order() {
            let g = &mut self.groups[gid as usize];
            if g.status == ReadGroupStatus::Duplicate {
                continue;
            }
            if self.expected_read_ids.remove(&g.id) {
                g.has_remote_junction_reads = true;
                g.status = ReadGroupStatus::Expected;
                g.expected = true;
                self.expected_groups.push(gid);
            }
        }
    }

    pub fn assign_junction_fragments_and_support(&mut self) {
        let permitted_position_diff = if self.cfg.filter.sequencing.is_sbx() { SBX_MAX_DUPLICATE_DISTANCE } else { 0 };
        self.duplicate_groups = self.mark_supplementary_duplicates(permitted_position_diff);

        let order = self.map_order();
        let mut candidate_support_groups: Vec<GroupId> = Vec::new();
        let has_expected = !self.expected_read_ids.is_empty();
        for &gid in &order {
            if self.groups[gid as usize].status == ReadGroupStatus::Duplicate {
                continue;
            }
            if has_expected {
                let g = &mut self.groups[gid as usize];
                if self.expected_read_ids.remove(&g.id) {
                    g.has_remote_junction_reads = true;
                    g.status = ReadGroupStatus::Expected;
                    g.expected = true;
                    self.expected_groups.push(gid);
                }
            }
            let g = &self.groups[gid as usize];
            // a short overlapping fragment is likely adapter
            if g.reads.iter().any(|&i| self.reads[i as usize].has_filter(F_INSERT_MAP_OVERLAP)) {
                continue;
            }
            if group_all_no_support(g, &self.reads) {
                continue;
            }
            if group_has_read_type(g, &self.reads, ReadType::Junction) {
                self.create_junction(gid);
                candidate_support_groups.push(gid);
            } else if group_has_read_type(g, &self.reads, ReadType::CandidateSupport) {
                candidate_support_groups.push(gid);
            }
        }
        for k in 0..self.junctions.len() {
            self.set_initial_read(self.junctions[k]);
        }
        self.recent_junction_index = INVALID_JUNC_INDEX;
        self.initial_supporting_frags = candidate_support_groups.len() as i32;

        let has_junctions = !self.junctions.is_empty();
        if has_junctions {
            // by the first read's start, so that the last junction index is a good place to look from
            candidate_support_groups.sort_by_key(|&g| self.reads[self.groups[g as usize].reads[0] as usize].start);
        }
        let has_paired_reads = self.cfg.filter.sequencing.has_paired_reads();
        let mut supported: Vec<(JuncId, ReadType)> = Vec::new();
        for &gid in &candidate_support_groups {
            supported.clear();
            if has_junctions {
                let mut k = 0;
                while k < self.groups[gid as usize].reads.len() {
                    let rid = self.groups[gid as usize].reads[k];
                    let t = self.reads[rid as usize].read_type;
                    if t == ReadType::CandidateSupport || t == ReadType::Junction {
                        self.check_junction_support(gid, rid, &mut supported);
                    }
                    k += 1;
                }
            }
            if !supported.is_empty() {
                for &(jid, t) in &supported {
                    let (pos, orient) = (self.jstore[jid as usize].position, self.jstore[jid as usize].orient);
                    if self.groups[gid as usize].has_junction_position(pos, orient) {
                        continue;
                    }
                    if t == ReadType::ExactSupport {
                        self.jstore[jid as usize].exact_support_groups.push(gid);
                    } else {
                        self.jstore[jid as usize].supporting_groups.push(gid);
                    }
                    self.groups[gid as usize].add_junction_position(pos, orient);
                }
            } else if !self.groups[gid as usize].has_junction_positions() {
                self.add_remote_candidate(gid);
            }
            let mut is_discordant_group = false;
            if has_paired_reads {
                let g = &self.groups[gid as usize];
                is_discordant_group = self.discordant_finder.is_discordant_group(g, &self.reads);
                if is_discordant_group && self.discordant_finder.is_relevant_discordant_group(g, &self.reads) {
                    self.candidate_discordant_groups.push(gid);
                }
            }
            let g = &self.groups[gid as usize];
            if is_discordant_group || DiscordantStats::has_primary_with_supplementary(g, &self.reads) {
                self.discordant_stats.process_read_group(g, &self.reads, is_discordant_group, &self.cfg.contigs);
            }
        }
    }

    fn add_remote_candidate(&mut self, gid: GroupId) {
        let g = &mut self.groups[gid as usize];
        if !g.remote_candidate {
            g.remote_candidate = true;
            self.remote_candidate_groups.push(gid);
        }
    }

    /// JunctionUtils.markSupplementaryDuplicates: among reads with a supplementary alignment that share an
    /// unclipped 5' position, a primary and a supplementary of fragments with the same other end are
    /// duplicates of each other and one of the two groups is marked.
    pub fn mark_supplementary_duplicates(&mut self, permitted_position_diff: i32) -> i32 {
        let is_candidate = |t: ReadType| t == ReadType::Junction || t == ReadType::CandidateSupport;
        let mut initial_position_map: JavaHashMap<i32, Vec<ReadId>> = JavaHashMap::new();
        for gid in self.group_map.values() {
            let g = &self.groups[*gid as usize];
            if !g.reads.iter().any(|&i| self.reads[i as usize].has_supp_alignment()) {
                continue;
            }
            if !g.reads.iter().any(|&i| is_candidate(self.reads[i as usize].read_type)) {
                continue;
            }
            for &rid in &g.reads {
                let read = &self.reads[rid as usize];
                if !read.has_supp_alignment() {
                    continue;
                }
                if read.is_paired() && !read.is_mate_mapped() {
                    continue;
                }
                if !is_candidate(read.read_type) {
                    continue;
                }
                let unclipped_position = if read.orientation().is_forward() { read.unclipped_start } else { read.unclipped_end };
                let key = if permitted_position_diff == 0 {
                    initial_position_map.contains_key(&unclipped_position).then_some(unclipped_position)
                } else {
                    initial_position_map.keys().copied().find(|k| (k - unclipped_position).abs() <= permitted_position_diff)
                };
                match key {
                    Some(k) => initial_position_map.get_mut(&k).unwrap().push(rid),
                    None => {
                        initial_position_map.insert(unclipped_position, vec![rid]);
                    }
                }
                break;
            }
        }

        let mut duplicate_fragment_count = 0;
        let mut duplicates: Vec<GroupId> = Vec::new();
        'lists: for list in initial_position_map.values() {
            if list.len() < 2 {
                continue;
            }
            let fragments: Vec<FragmentData> = list.iter().map(|&i| FragmentData::new(i, &self.reads[i as usize], &self.cfg.contigs)).collect();
            if !fragments.iter().any(|f| !f.is_primary) || !fragments.iter().any(|f| f.is_primary) {
                continue;
            }
            for i in 0..fragments.len() - 1 {
                let first = &fragments[i];
                for next in &fragments[i + 1..] {
                    if first.is_primary == next.is_primary {
                        continue;
                    }
                    if permitted_position_diff == 0 {
                        if !first.matches(next) {
                            continue;
                        }
                    } else if !first.within_position_range(next, permitted_position_diff) {
                        continue;
                    }
                    // the supplementary goes when it is the lower of the two, but a consensus read is kept
                    let mark_next_as_duplicate =
                        if first.is_consensus != next.is_consensus { first.is_consensus } else { first.is_primary == first.read_is_lower_vs_supp };
                    let dup_read = if mark_next_as_duplicate { next.read } else { first.read };
                    let dup_id: Box<[u8]> = self.trimmer.trim(self.reads[dup_read as usize].id()).into();
                    match self.group_map.get(&dup_id) {
                        Some(&gid) => {
                            duplicates.push(gid);
                            duplicate_fragment_count += 1;
                        }
                        None => {
                            crate::error!(
                                "duplicateReadId({}) group not found from read({}), trimmerEnabled({})",
                                String::from_utf8_lossy(&dup_id),
                                String::from_utf8_lossy(self.reads[dup_read as usize].id()),
                                self.trimmer.enabled()
                            );
                            break 'lists;
                        }
                    }
                }
            }
        }
        for gid in duplicates {
            self.groups[gid as usize].status = ReadGroupStatus::Duplicate;
        }
        duplicate_fragment_count
    }

    fn create_junction(&mut self, gid: GroupId) {
        let mut junctions: Vec<JuncId> = Vec::with_capacity(2);
        let mut k = 0;
        while k < self.groups[gid as usize].reads.len() {
            let rid = self.groups[gid as usize].reads[k];
            k += 1;
            let read = &self.reads[rid as usize];
            if read.is_unmapped() || read.read_type != ReadType::Junction {
                continue;
            }
            if let Some(indel) = read.indel_coords() {
                self.handle_indel_junction(gid, rid, indel);
            }
            let read = &self.reads[rid as usize];
            if read.has_filter(F_SOFT_CLIP_LENGTH) {
                continue;
            }
            for is_start in [true, false] {
                let read = &self.reads[rid as usize];
                let sc_length = if is_start { read.sc_left } else { read.sc_right };
                // the soft-clip filter has allowed for LINE reads, so the shorter LINE length applies here
                if sc_length < MIN_LINE_SOFT_CLIP_LENGTH {
                    continue;
                }
                let orientation = if is_start { Orientation::Reverse } else { Orientation::Forward };
                let position = if orientation.is_reverse() { read.start } else { read.end };
                if self.region.contains_position(position) {
                    let jid = self.get_or_create_junction(rid, position, orientation);
                    self.jstore[jid as usize].add_read_type(rid, ReadType::Junction);
                    if !junctions.contains(&jid) {
                        junctions.push(jid);
                    }
                }
            }
        }
        for &jid in &junctions {
            self.jstore[jid as usize].junction_groups.push(gid);
        }
        for &jid in &junctions {
            let (p, o) = (self.jstore[jid as usize].position, self.jstore[jid as usize].orient);
            self.groups[gid as usize].add_junction_position(p, o);
        }
    }

    fn handle_indel_junction(&mut self, gid: GroupId, rid: ReadId, indel: IndelCoords) {
        if indel.length < self.cfg.filter.min_indel_length {
            return;
        }
        // enough aligned bases on both sides to anchor the indel
        {
            let cigar = self.reads[rid as usize].cigar();
            let last_index = cigar.len() as i32 - 1;
            for (i, &c) in cigar.iter().enumerate() {
                let i = i as i32;
                let short_match = op(c) == OP_M && oplen(c) < MIN_INDEL_ANCHOR_LENGTH;
                if (i <= 1 || i >= last_index - 1) && short_match {
                    return;
                }
            }
        }
        let junction_start = self.get_or_create_junction(rid, indel.pos_start, Orientation::Forward);
        let junction_end = self.get_or_create_junction(rid, indel.pos_end, Orientation::Reverse);
        self.jstore[junction_start as usize].linked_indel = Some(junction_end);
        self.jstore[junction_end as usize].linked_indel = Some(junction_start);
        for jid in [junction_start, junction_end] {
            let j = &mut self.jstore[jid as usize];
            j.internal_indel = true;
            j.junction_groups.push(gid);
            j.add_read_type(rid, ReadType::Junction);
            let (p, o) = (j.position, j.orient);
            self.groups[gid as usize].add_junction_position(p, o);
        }
    }

    fn find_index(&self, position: i32) -> i32 {
        find_junction_index(self.junctions.len(), |i| self.jstore[self.junctions[i] as usize].position, position, SIMPLE_SEARCH_COUNT)
    }

    /// JunctionUtils.findJunctionMatchIndex: the same position and orientation at the index or next to it
    fn find_match_index(&self, position: i32, orientation: Orientation, junc_index: i32) -> i32 {
        for i in -1..=1 {
            let index = junc_index + i;
            if index >= 0 && (index as usize) < self.junctions.len() {
                let j = &self.jstore[self.junctions[index as usize] as usize];
                if j.position == position && j.orient == orientation {
                    return index;
                }
            }
        }
        INVALID_JUNC_INDEX
    }

    fn get_or_create_junction(&mut self, rid: ReadId, junction_position: i32, orientation: Orientation) -> JuncId {
        if self.recent_junction_index >= 0 {
            let jid = self.junctions[self.recent_junction_index as usize];
            let j = &self.jstore[jid as usize];
            if j.position == junction_position && j.orient == orientation {
                return jid;
            }
        }
        let close = self.find_index(junction_position);
        if close != INVALID_JUNC_INDEX {
            let exact = self.find_match_index(junction_position, orientation, close);
            if exact != INVALID_JUNC_INDEX {
                self.recent_junction_index = exact;
                return self.junctions[exact as usize];
            }
        }
        let new_index = close + 1;
        let jid = self.jstore.len() as JuncId;
        self.jstore.push(JunctionData::new(junction_position, orientation, Some(rid)));
        if new_index >= 0 && (new_index as usize) < self.junctions.len() {
            self.junctions.insert(new_index as usize, jid);
        } else {
            self.junctions.push(jid);
        }
        self.recent_junction_index = new_index;
        jid
    }

    /// addJunctionData of the Java tests: a junction at the read's clipped side
    #[cfg(test)]
    pub fn add_junction_data(&mut self, rid: ReadId) -> JuncId {
        let read = &self.reads[rid as usize];
        let orientation = if read.is_left_clipped() { Orientation::Reverse } else { Orientation::Forward };
        let position = if orientation.is_reverse() { read.start } else { read.end };
        self.get_or_create_junction(rid, position, orientation)
    }

    /// Adds a read to the store without a group, for tests that build junctions directly.
    pub fn add_read(&mut self, read: PrepRead) -> ReadId {
        self.reads.push(read);
        (self.reads.len() - 1) as ReadId
    }

    /// JunctionData.setInitialRead: the junction read with the most high-quality soft-clip bases
    fn set_initial_read(&mut self, jid: JuncId) {
        let j = &self.jstore[jid as usize];
        if j.internal_indel {
            return;
        }
        if j.junction_reads.len() == 1 {
            let r = j.junction_reads[0];
            self.jstore[jid as usize].top_junction_read = Some(r);
            return;
        }
        let use_left = j.orient.is_reverse();
        let mut max_high_qual_bases = 0i32;
        let mut top: Option<ReadId> = None;
        for &rid in &j.junction_reads {
            let read = &self.reads[rid as usize];
            let sc_length = if use_left { read.sc_left } else { read.sc_right };
            if sc_length < max_high_qual_bases {
                continue;
            }
            let quals = read.base_qualities();
            let n = quals.len() as i32;
            let (range_start, range_end) = if use_left { (0, sc_length) } else { (n - sc_length, n) };
            if range_end - range_start + 1 < max_high_qual_bases {
                continue;
            }
            let t = crate::common::base_qual::low_base_qual_threshold();
            let above = quals[range_start as usize..range_end as usize].iter().filter(|&&q| !below_qual_threshold(q, t)).count() as i32;
            if above > max_high_qual_bases {
                top = Some(rid);
                max_high_qual_bases = above;
            }
        }
        if top.is_some() {
            self.jstore[jid as usize].top_junction_read = top;
        }
    }

    fn supported_get(supported: &[(JuncId, ReadType)], jid: JuncId) -> Option<ReadType> {
        supported.iter().find(|x| x.0 == jid).map(|x| x.1)
    }

    fn supported_put(supported: &mut Vec<(JuncId, ReadType)>, jid: JuncId, t: ReadType) {
        match supported.iter_mut().find(|x| x.0 == jid) {
            Some(x) => x.1 = t,
            None => supported.push((jid, t)),
        }
    }

    fn check_junction_support(&mut self, gid: GroupId, rid: ReadId, supported: &mut Vec<(JuncId, ReadType)>) {
        let has_paired_reads = self.cfg.filter.sequencing.has_paired_reads();
        let distance_buffer = if has_paired_reads { self.cfg.filter.max_supporting_fragment_distance } else { UNPAIRED_READ_JUNCTION_DISTANCE };
        let (read_start, read_end, unclipped_start, unclipped_end) = {
            let r = &self.reads[rid as usize];
            (r.start, r.end, r.unclipped_start, r.unclipped_end)
        };
        let within = |pos: i32| (read_end - pos).abs() <= distance_buffer || (read_start - pos).abs() <= distance_buffer;

        // the last index first, since the next read is likely to be close by
        let mut close = INVALID_JUNC_INDEX;
        if self.recent_junction_index >= 0 {
            let recent = self.recent_junction_index as usize;
            let pos = self.jstore[self.junctions[recent] as usize].position;
            if within(pos) || position_within(pos, read_start, read_end) {
                close = self.recent_junction_index;
            } else if recent < self.junctions.len() - 1 {
                let next_pos = self.jstore[self.junctions[recent + 1] as usize].position;
                if positions_overlap(read_start, read_end, pos, next_pos) {
                    close = self.recent_junction_index;
                }
            }
        }
        if close == INVALID_JUNC_INDEX {
            close = self.find_index(read_start);
        }
        if close == INVALID_JUNC_INDEX {
            close = 0;
        }
        self.recent_junction_index = close;

        for search_up in [true, false] {
            let mut index = if search_up { close } else { close - 1 };
            while index >= 0 && (index as usize) < self.junctions.len() {
                let jid = self.junctions[index as usize];
                let pos = self.jstore[jid as usize].position;
                if !has_paired_reads {
                    if search_up && pos >= unclipped_end {
                        break;
                    }
                    if !search_up && pos < unclipped_start {
                        break;
                    }
                    if !within(pos) {
                        index += if search_up { 1 } else { -1 };
                        continue;
                    }
                } else if !within(pos) {
                    // the find routine gives a match or the junction before, so look at least one higher
                    if !search_up || index > close {
                        break;
                    }
                }
                self.check_read_supports_junction(gid, rid, jid, supported);
                index += if search_up { 1 } else { -1 };
            }
        }
        self.check_indel_support(rid, supported, close);
    }

    fn check_read_supports_junction(&mut self, gid: GroupId, rid: ReadId, jid: JuncId, supported: &mut Vec<(JuncId, ReadType)>) {
        let j = &self.jstore[jid as usize];
        if self.groups[gid as usize].has_junction_position(j.position, j.orient) {
            return;
        }
        let read_type = Self::supported_get(supported, jid);
        let read = &self.reads[rid as usize];
        if read_type != Some(ReadType::ExactSupport) && has_exact_junction_support(read, j, j.top_junction_read.map(|t| &self.reads[t as usize])) {
            self.jstore[jid as usize].add_read_type(rid, ReadType::ExactSupport);
            self.reads[rid as usize].set_read_type_ranked(ReadType::ExactSupport);
            Self::supported_put(supported, jid, ReadType::ExactSupport);
            return;
        }
        if read_type != Some(ReadType::Support) && self.cfg.filter.sequencing.has_paired_reads() && has_other_junction_support(read, j, &self.cfg.filter) {
            self.jstore[jid as usize].add_read_type(rid, ReadType::Support);
            self.reads[rid as usize].set_read_type_ranked(ReadType::Support);
            Self::supported_put(supported, jid, ReadType::Support);
        }
    }

    fn check_indel_support(&mut self, rid: ReadId, supported: &mut Vec<(JuncId, ReadType)>, close: i32) {
        let read = &self.reads[rid as usize];
        let Some(indel) = read.indel_coords() else { return };
        let (mut implied_start, mut implied_end) = (read.start, read.end);
        if indel.is_insert() {
            implied_start -= indel.length;
            implied_end += indel.length;
        } else {
            implied_start += indel.length;
            implied_end -= indel.length;
        }
        let read_bounds_min = read.start.min(implied_start);
        let read_bounds_max = read.end.max(implied_end);
        // a read with a long enough indel need only cover a junction with any of its bases
        for search_up in [true, false] {
            let mut index = if search_up { close } else { close - 1 };
            while index >= 0 && (index as usize) < self.junctions.len() {
                let jid = self.junctions[index as usize];
                let (pos, is_reverse) = (self.jstore[jid as usize].position, self.jstore[jid as usize].is_reverse());
                if search_up && pos > read_bounds_max {
                    break;
                }
                if !search_up && read_bounds_min > pos {
                    break;
                }
                if Self::supported_get(supported, jid).is_some() {
                    index += if search_up { 1 } else { -1 };
                    continue;
                }
                let mut has_exact_support = false;
                for is_start in [true, false] {
                    let indel_pos = if is_start { indel.pos_start } else { indel.pos_end };
                    if indel_pos != pos {
                        continue;
                    }
                    if is_start == is_reverse {
                        continue;
                    }
                    self.reads[rid as usize].set_read_type_ranked(ReadType::ExactSupport);
                    self.jstore[jid as usize].add_read_type(rid, ReadType::ExactSupport);
                    Self::supported_put(supported, jid, ReadType::ExactSupport);
                    has_exact_support = true;
                }
                if !has_exact_support {
                    self.reads[rid as usize].set_read_type_ranked(ReadType::Support);
                    self.jstore[jid as usize].add_read_type(rid, ReadType::Support);
                    Self::supported_put(supported, jid, ReadType::Support);
                }
                index += if search_up { 1 } else { -1 };
            }
        }
    }

    fn filter_junctions(&mut self, depth: &DepthTracker) {
        let mut removed_groups: Vec<GroupId> = Vec::new();
        let mut index = 0;
        while index < self.junctions.len() {
            let jid = self.junctions[index];
            self.jstore[jid as usize].depth = depth.calc_depth(self.jstore[jid as usize].position) as i32;
            let mut has_support = self.junction_has_support(jid);
            if !has_support {
                if let Some(linked) = self.jstore[jid as usize].linked_indel {
                    has_support = self.junction_has_support(linked);
                }
            }
            if has_support {
                index += 1;
            } else {
                self.junctions.remove(index);
                let j = &self.jstore[jid as usize];
                removed_groups.extend(j.junction_groups.iter().chain(&j.exact_support_groups).chain(&j.supporting_groups));
            }
        }
        // the junction positions of the groups, without those of the junctions just dropped
        for &gid in self.group_map.values() {
            self.groups[gid as usize].junction_positions.clear();
        }
        for &jid in &self.junctions {
            let j = &self.jstore[jid as usize];
            for &gid in j.junction_groups.iter().chain(&j.exact_support_groups).chain(&j.supporting_groups) {
                self.groups[gid as usize].add_junction_position(j.position, j.orient);
            }
        }
        // groups left without a junction go back to being candidates, for the spanning cache
        let mut seen = vec![false; self.groups.len()];
        for gid in removed_groups {
            if std::mem::replace(&mut seen[gid as usize], true) {
                continue;
            }
            if self.groups[gid as usize].status == ReadGroupStatus::Duplicate {
                for &rid in &self.groups[gid as usize].reads {
                    self.reads[rid as usize].read_type = ReadType::NoSupport;
                }
                continue;
            }
            if !self.groups[gid as usize].has_junction_positions() {
                self.add_remote_candidate(gid);
                for &rid in &self.groups[gid as usize].reads {
                    let r = &mut self.reads[rid as usize];
                    if r.read_type != ReadType::NoSupport {
                        r.read_type = ReadType::CandidateSupport;
                    }
                }
            }
        }
    }

    fn junction_has_support(&mut self, jid: JuncId) -> bool {
        let j = &self.jstore[jid as usize];
        if !j.hotspot && !self.junction_above_min_af(j) {
            return false;
        }
        if j.discordant_group {
            return true;
        }
        let junction_frags = j.junction_fragment_count();
        let exact_support_count = j.exact_support_fragment_count();
        if self.junction_matches_hotspot(j) {
            self.jstore[jid as usize].hotspot = true;
            if junction_frags + exact_support_count >= MIN_HOTSPOT_JUNCTION_SUPPORT {
                return true;
            }
        }
        if !self.jstore[jid as usize].internal_indel {
            match self.find_well_anchored_read_soft_clip_length(jid) {
                None => return false,
                Some(len) => self.jstore[jid as usize].max_valid_soft_clip_length = len,
            }
        }
        junction_frags + exact_support_count >= self.cfg.filter.min_junction_support
    }

    /// Sites with support below 0.5% of the depth (1% for discordant groups, 2% for other sequencing) are dropped.
    fn junction_above_min_af(&self, j: &JunctionData) -> bool {
        const DEPTH_MIN_CHECK: i32 = 100;
        let region_depth = j.depth;
        if region_depth < DEPTH_MIN_CHECK {
            return true;
        }
        let required_ratio = if j.discordant_group {
            0.01
        } else if self.cfg.filter.sequencing.is_illumina() {
            0.005
        } else {
            0.02
        };
        let mut junction_support = j.junction_fragment_count() + j.exact_support_fragment_count();
        if j.discordant_group {
            junction_support += j.supporting_fragment_count();
        }
        junction_support as f64 >= region_depth as f64 * required_ratio
    }

    /// KnownHotspot.junctionMatchesHotspot: a junction in one side of a hotspot with a mate in the other
    fn junction_matches_hotspot(&self, j: &JunctionData) -> bool {
        let group_matches_region = |gid: GroupId, region: &Region| {
            self.groups[gid as usize].reads.iter().any(|&rid| {
                let r = &self.reads[rid as usize];
                region.contains(r.mate_chr, mate_start(&r.rec))
            })
        };
        for h in &self.known_hotspots {
            if h.region_start.contains_position(j.position) {
                if j.junction_groups.iter().any(|&g| group_matches_region(g, &h.region_end)) {
                    return true;
                }
            } else if h.region_end.contains_position(j.position) && j.junction_groups.iter().any(|&g| group_matches_region(g, &h.region_start)) {
                return true;
            }
        }
        false
    }

    /// JunctionUtils.findWellAnchoredReadSoftClipLength: the longest soft clip of a junction read whose
    /// alignment holds up once repeats are taken off. Sorts the junction's reads by soft-clip length.
    fn find_well_anchored_read_soft_clip_length(&mut self, jid: JuncId) -> Option<i32> {
        let forward = self.jstore[jid as usize].is_forward();
        let mut junction_reads = std::mem::take(&mut self.jstore[jid as usize].junction_reads);
        junction_reads.sort_by_key(|&r| if forward { -self.reads[r as usize].sc_right } else { -self.reads[r as usize].sc_left });
        let mut result = None;
        for &rid in &junction_reads {
            let read = &self.reads[rid as usize];
            let valid = if read.has_line_tail {
                true
            } else {
                let score = calc_repeat_trimmed_alignment_score_read(read, MIN_CALC_ALIGNMENT_SCORE, true);
                if score < 0.0 {
                    continue;
                }
                score >= MIN_CALC_ALIGNMENT_SCORE as f64
                    || (score >= MIN_CALC_ALIGNMENT_LOWER_SCORE as f64 && above_aligned_score_difference(read, self.cfg.filter.min_alignment_bases))
            };
            if valid {
                result = Some(if forward { read.sc_right } else { read.sc_left });
                break;
            }
        }
        self.jstore[jid as usize].junction_reads = junction_reads;
        result
    }

    pub fn find_discordant_groups(&mut self) {
        if !self.cfg.filter.sequencing.has_paired_reads() {
            return;
        }
        let candidates = std::mem::take(&mut self.candidate_discordant_groups);
        let first = self.jstore.len();
        self.discordant_finder.form_discordant_junctions(&candidates, &mut self.groups, &mut self.reads, &mut self.jstore);
        if candidates.len() > 2000 && self.jstore.len() > first {
            crate::debug!(
                "region({}) found {} discordant group junctions from {} read groups",
                self.region.display(&self.cfg.contigs),
                self.jstore.len() - first,
                candidates.len()
            );
        }
        self.candidate_discordant_groups = candidates;
        let (outside, other_chr) = (self.discordant_finder.remote_outside_partition.get(), self.discordant_finder.remote_other_chromosome.get());
        if outside > 0 {
            crate::debug!(
                "region({}) {} discordant junctions are at remote regions outside the partition, {} on another chromosome",
                self.region.display(&self.cfg.contigs),
                outside,
                other_chr
            );
        }
        let java_p8 = crate::common::compat::java("P8");
        for jid in first..self.jstore.len() {
            let j = &self.jstore[jid];
            if let (false, Some(chr)) = (java_p8, j.outside_chr) {
                // Java: the junction joins this partition's list, is matched against the junctions here by
                // position alone and is written under this partition's chromosome
                let (position, orient) = (j.position, j.orient);
                let known = self.outside_junctions.iter().any(|&o| {
                    let x = &self.jstore[o as usize];
                    x.outside_chr == Some(chr) && x.position == position && x.orient == orient
                });
                if !known {
                    self.outside_junctions.push(jid as JuncId);
                }
            } else {
                self.add_discordant_junction(jid as JuncId);
            }
        }
    }

    /// addDiscordantJunction: into the ordered list, in place of a junction at the same site that is neither a
    /// hotspot nor has the fragments to stand by itself
    pub fn add_discordant_junction(&mut self, jid: JuncId) {
        let (position, orient, is_disc) = {
            let j = &self.jstore[jid as usize];
            (j.position, j.orient, j.discordant_group)
        };
        let close = self.find_index(position);
        if close != INVALID_JUNC_INDEX {
            let exact = self.find_match_index(position, orient, close);
            if exact != INVALID_JUNC_INDEX {
                let existing = &self.jstore[self.junctions[exact as usize] as usize];
                if !existing.discordant_group && is_disc && !existing.hotspot && existing.junction_fragment_count() < self.cfg.filter.min_junction_support {
                    self.junctions[exact as usize] = jid;
                }
                return;
            }
        }
        self.junctions.insert((close + 1) as usize, jid);
    }

    fn match_junctions_to_saga(&mut self) {
        let Some(saga) = &self.saga else { return };
        let chr = self.cfg.contigs.name(self.region.chr).to_string();
        for &jid in &self.junctions {
            let j = &mut self.jstore[jid as usize];
            j.saga_match = saga.match_junction(&chr, j.position, j.orient);
        }
        // the matcher covers this chromosome around the partition only
        for &jid in &self.outside_junctions {
            let j = &mut self.jstore[jid as usize];
            if j.outside_chr == Some(self.region.chr) {
                j.saga_match = saga.match_junction(&chr, j.position, j.orient);
            }
        }
    }

    /// formUniqueAssignedGroups: the groups of the junctions, each once, then the expected groups
    pub fn form_unique_assigned_groups(&mut self) -> Vec<GroupId> {
        let mut junction_groups: Vec<GroupId> = Vec::new();
        // Java keeps a set of read ids; the groups of the map have distinct ids
        let mut seen = vec![false; self.groups.len()];
        for &jid in self.junctions.iter().chain(&self.outside_junctions) {
            let j = &self.jstore[jid as usize];
            for &gid in j.junction_groups.iter().chain(&j.supporting_groups).chain(&j.exact_support_groups) {
                if !std::mem::replace(&mut seen[gid as usize], true) {
                    junction_groups.push(gid);
                }
            }
        }
        for k in 0..self.expected_groups.len() {
            let gid = self.expected_groups[k];
            if seen[gid as usize] {
                continue;
            }
            self.groups[gid as usize].status = ReadGroupStatus::Expected;
            for &rid in &self.groups[gid as usize].reads {
                self.reads[rid as usize].read_type = ReadType::Expected;
            }
            junction_groups.push(gid);
        }
        junction_groups
    }

    /// getRemoteCandidateReadGroups: groups with a read in another partition and no junction, other than the
    /// expected ones
    pub fn get_remote_candidate_read_groups(&self) -> Vec<GroupId> {
        self.remote_candidate_groups
            .iter()
            .copied()
            .filter(|&g| !self.groups[g as usize].expected)
            .filter(|&g| self.groups[g as usize].junction_positions.is_empty())
            .collect()
    }
}
