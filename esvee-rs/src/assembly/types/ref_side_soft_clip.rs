//! RefSideSoftClip: reads that soft-clip on the reference side of a junction, suggesting a short templated insertion.

use crate::assembly::constants::*;
use crate::assembly::orientation::Orientation;
use crate::assembly::ptr::P;
use crate::assembly::read::read::Read;
use crate::assembly::types::junction::Junction;
use crate::assembly::types::junction_assembly::AsmP;
use crate::common::java_hash::JavaHashSet;
use crate::common::java_sort::java_sort;

pub type RefSideSoftClipP = P<RefSideSoftClip>;

pub struct RefSideSoftClip {
    pub position: i32,
    pub orient: Orientation,
    read_ids: JavaHashSet<String>,
    max_length: i32,
    matched_original: bool,
}

impl RefSideSoftClip {
    pub fn new(position: i32, orientation: Orientation, read: &Read, read_soft_clip_length: i32) -> RefSideSoftClip {
        let mut r = RefSideSoftClip { position, orient: orientation, read_ids: JavaHashSet::new(), max_length: 0, matched_original: false };
        r.add_read(read, read_soft_clip_length);
        r
    }
    pub fn read_ids(&self) -> &JavaHashSet<String> {
        &self.read_ids
    }
    pub fn max_length(&self) -> i32 {
        self.max_length
    }
    pub fn read_count(&self) -> i32 {
        self.read_ids.len() as i32
    }
    pub fn mark_matches_original(&mut self) {
        self.matched_original = true;
    }
    pub fn matches_original(&self) -> bool {
        self.matched_original
    }
    pub fn add_read(&mut self, read: &Read, read_soft_clip_length: i32) {
        self.read_ids.insert(read.id().to_string());
        self.max_length = self.max_length.max(read_soft_clip_length);
    }
    pub fn has_proximate_match(&self, other_ref_position: i32) -> bool {
        (self.position - other_ref_position).abs() <= PROXIMATE_REF_SIDE_SOFT_CLIPS
    }

    pub fn check_add_ref_side_soft_clip(ref_side_soft_clips: &mut Vec<RefSideSoftClipP>, junction: &Junction, read: &Read) -> bool {
        let ref_side_soft_clip_position;
        let ref_side_soft_clip_length;
        if junction.is_forward() {
            ref_side_soft_clip_position = read.alignment_start();
            if ref_side_soft_clip_position > junction.position - PHASED_ASSEMBLY_MIN_TI {
                return false;
            }
            ref_side_soft_clip_length = read.left_clip_length();
        } else {
            ref_side_soft_clip_position = read.alignment_end();
            if ref_side_soft_clip_position < junction.position + PHASED_ASSEMBLY_MIN_TI {
                return false;
            }
            ref_side_soft_clip_length = read.right_clip_length();
        }
        if ref_side_soft_clip_length < REF_SIDE_INITIAL_MIN_SOFT_CLIP_LENGTH {
            return false;
        }
        match ref_side_soft_clips.iter().find(|x| x.position == ref_side_soft_clip_position) {
            None => ref_side_soft_clips.push(P::new(RefSideSoftClip::new(ref_side_soft_clip_position, junction.orient.opposite(), read, ref_side_soft_clip_length))),
            Some(existing) => existing.m().add_read(read, ref_side_soft_clip_length),
        }
        true
    }

    pub fn purge_ref_side_soft_clips(ref_side_soft_clips: &mut Vec<RefSideSoftClipP>, non_soft_clip_ref_position: i32) {
        Self::purge_ref_side_soft_clips_ext(ref_side_soft_clips, ASSEMBLY_MIN_READ_SUPPORT, REF_SIDE_MIN_SOFT_CLIP_LENGTH, non_soft_clip_ref_position)
    }

    pub fn purge_ref_side_soft_clips_ext(ref_side_soft_clips: &mut Vec<RefSideSoftClipP>, min_count: i32, min_length: i32, non_soft_clip_ref_position: i32) {
        if ref_side_soft_clips.is_empty() {
            return;
        }
        let mut index = 0;
        let mut matching: Option<RefSideSoftClipP> = None;
        while index < ref_side_soft_clips.len() {
            let sc = ref_side_soft_clips[index];
            let is_position_match_or_close = sc.has_proximate_match(non_soft_clip_ref_position);
            let insufficient_support_or_length = sc.read_count() < min_count || sc.max_length() < min_length;
            let invalid_position = if sc.orient.is_forward() { sc.position > non_soft_clip_ref_position } else { sc.position < non_soft_clip_ref_position };
            if insufficient_support_or_length || is_position_match_or_close || invalid_position {
                ref_side_soft_clips.remove(index);
                if is_position_match_or_close {
                    if sc.position == non_soft_clip_ref_position {
                        matching = Some(sc);
                    } else if matching.is_none() || matching.unwrap().read_count() < sc.read_count() {
                        matching = Some(sc);
                    }
                }
            } else {
                index += 1;
            }
        }
        if let Some(m) = matching {
            m.m().mark_matches_original();
            ref_side_soft_clips.push(m);
        }
    }

    pub fn check_support_vs_ref_side_soft_clip(assembly: AsmP) {
        let mut ref_side_soft_clips: Vec<RefSideSoftClipP> = Vec::new();
        for read in assembly.support().iter() {
            Self::check_add_ref_side_soft_clip(&mut ref_side_soft_clips, &assembly.junction(), read.cached_read().expect("support read without cached read"));
        }
        if ref_side_soft_clips.is_empty() {
            return;
        }
        java_sort(&mut ref_side_soft_clips, |a, b| (-a.read_count()).cmp(&(-b.read_count())));
        let main = ref_side_soft_clips[0];
        if main.read_count() < ASSEMBLY_MIN_READ_SUPPORT {
            return;
        }
        let main_soft_clip_count = main.read_count();
        let total_soft_clip_count: i32 = ref_side_soft_clips.iter().map(|x| x.read_count()).sum();
        let non_soft_clip_count = assembly.support_count() as i32 - total_soft_clip_count;
        if non_soft_clip_count > 0 && non_soft_clip_count < main_soft_clip_count && non_soft_clip_count < ASSEMBLY_SPLIT_MIN_READ_SUPPORT {
            let support = assembly.m().support_mut();
            let mut index = 0;
            while index < support.len() {
                let keep = ref_side_soft_clips.iter().any(|x| x.read_ids().contains(&support[index].id().to_string()));
                if !keep {
                    support.remove(index);
                } else {
                    index += 1;
                }
            }
        }
    }
}

impl std::fmt::Display for RefSideSoftClip {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} len({}) reads({}) {}", self.position, self.max_length, self.read_ids.len(), if self.matched_original { "matchOrig" } else { "" })
    }
}
