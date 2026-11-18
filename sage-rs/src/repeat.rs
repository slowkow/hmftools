//! RepeatInfo and RepeatBoundaries.

use crate::constants::*;
use crate::types::{position_within, positions_overlap};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RepeatInfo {
    pub index: i32,
    pub bases: Vec<u8>,
    pub count: i32,
}

impl RepeatInfo {
    pub fn new(index: i32, bases: Vec<u8>, count: i32) -> Self {
        RepeatInfo { index, bases, count }
    }
    pub fn repeat_length(&self) -> i32 {
        self.bases.len() as i32
    }
    pub fn total_length(&self) -> i32 {
        self.count * self.repeat_length()
    }
    pub fn end_index(&self) -> i32 {
        self.index + self.total_length() - 1
    }
    pub fn bases_str(&self) -> &str {
        std::str::from_utf8(&self.bases).unwrap()
    }
    fn matches(&self, o: &RepeatInfo) -> bool {
        self.index == o.index && self.bases == o.bases && self.count == o.count
    }
    fn is_multiple_of(&self, other: &RepeatInfo) -> bool {
        if other.repeat_length() < self.repeat_length() {
            return false;
        }
        if other.index < self.index || other.end_index() > self.end_index() {
            return false;
        }
        let mut j = 0;
        for i in 0..other.bases.len() {
            if j >= self.bases.len() {
                j = 0;
            }
            if other.bases[i] != self.bases[j] {
                return false;
            }
            j += 1;
        }
        true
    }
    fn is_alternative_start(&self, other: &RepeatInfo) -> bool {
        let rl = self.repeat_length();
        if rl == 1 || other.repeat_length() != rl {
            return false;
        }
        if other.total_length() > self.total_length() {
            return false;
        }
        if !position_within(other.index, self.index, self.index + rl - 1) {
            return false;
        }
        let rl = rl as usize;
        for s in 1..rl {
            let mut all = true;
            let mut j = s;
            for i in 0..rl {
                if j >= rl {
                    j = 0;
                }
                if other.bases[i] != self.bases[j] {
                    all = false;
                    break;
                }
                j += 1;
            }
            if all {
                return true;
            }
        }
        false
    }
}

/// RepeatInfo.findMultiBaseRepeat; `repeat_len` is the unit length.
pub fn find_multi_base_repeat(bases: &[u8], index: i32, repeat_len: i32, min_count: i32) -> Option<RepeatInfo> {
    let n = bases.len() as i32;
    if index + (min_count * repeat_len) - 1 >= n {
        return None;
    }
    let mut i = index + repeat_len;
    while i < n - (repeat_len - 1) {
        let mut matched = 0;
        for j in 0..repeat_len {
            if bases[(i + j) as usize] == bases[(index + j) as usize] {
                matched += 1;
            } else {
                break;
            }
        }
        if matched != repeat_len {
            break;
        }
        i += repeat_len;
    }
    let count = (i - index) / repeat_len;
    if count < min_count {
        return None;
    }
    Some(RepeatInfo::new(index, bases[index as usize..(index + repeat_len) as usize].to_vec(), count))
}

pub fn extend_repeat_lower(repeat: RepeatInfo, bases: &[u8]) -> RepeatInfo {
    let mut extra = 0;
    let rl = repeat.repeat_length();
    let mut start = repeat.index - rl;
    while start > 0 {
        let mut matched = true;
        for i in 0..rl {
            if bases[(start + i) as usize] != repeat.bases[i as usize] {
                matched = false;
                break;
            }
        }
        if !matched {
            if extra > 0 {
                start += rl;
            }
            break;
        }
        extra += 1;
        if start - rl <= 0 {
            break;
        }
        start -= rl;
    }
    if extra == 0 {
        return repeat;
    }
    RepeatInfo::new(start, repeat.bases, repeat.count + extra)
}

fn add_if_unique(repeats: &mut Vec<RepeatInfo>, new: RepeatInfo) {
    for r in repeats.iter() {
        if new.matches(r) {
            return;
        }
        if new.is_multiple_of(r) || r.is_multiple_of(&new) {
            return;
        }
        if new.is_alternative_start(r) || r.is_alternative_start(&new) {
            return;
        }
    }
    repeats.push(new);
}

/// RepeatInfo.findMaxRepeat
pub fn find_max_repeat(
    bases: &[u8], search_start: i32, search_end: i32, max_length: i32, min_count: i32, extend_lower: bool,
    required_index: i32,
) -> Option<RepeatInfo> {
    let mut max_repeat: Option<RepeatInfo> = None;
    let mut index = search_start;
    while index <= (bases.len() as i32 - min_count).min(search_end) {
        for rl in 1..=max_length {
            let Some(mut repeat) = find_multi_base_repeat(bases, index, rl, min_count) else { continue };
            if extend_lower {
                repeat = extend_repeat_lower(repeat, bases);
            }
            if required_index >= 0 && !position_within(required_index, repeat.index, repeat.end_index()) {
                continue;
            }
            if max_repeat.as_ref().map_or(true, |m| repeat.count > m.count) {
                max_repeat = Some(repeat);
            }
        }
        index += 1;
    }
    max_repeat
}

pub struct RepeatBoundaries {
    pub lower_index: i32,
    pub upper_index: i32,
    pub max_repeat: Option<RepeatInfo>,
    pub all_repeats: Option<Vec<RepeatInfo>>,
}

const REPEAT_SEARCH_LENGTH: i32 = MAX_REPEAT_LENGTH * MIN_REPEAT_COUNT;

fn find_post_repeat_index(repeat: &RepeatInfo, bases: &[u8], search_down: bool) -> i32 {
    let mut partial = 0;
    let rl = repeat.repeat_length();
    if search_down {
        let mut base_index = repeat.index - 1;
        let mut in_repeat = rl - 1;
        while base_index >= 0 && partial < rl {
            if bases[base_index as usize] != repeat.bases[in_repeat as usize] {
                break;
            }
            partial += 1;
            base_index -= 1;
            in_repeat -= 1;
        }
        repeat.index - partial - 1
    } else {
        let mut base_index = repeat.end_index() + 1;
        let mut in_repeat = 0;
        while base_index < bases.len() as i32 && partial < rl {
            if bases[base_index as usize] != repeat.bases[in_repeat as usize] {
                break;
            }
            partial += 1;
            base_index += 1;
            in_repeat += 1;
        }
        repeat.end_index() + partial + 1
    }
}

fn find_extra_lower_repeat(bases: &[u8], required_start: i32) -> Option<RepeatInfo> {
    let mut index = 0.max(required_start - REPEAT_SEARCH_LENGTH);
    let mut max_repeat: Option<RepeatInfo> = None;
    while index <= (bases.len() as i32 - OUTER_MIN_REPEAT_COUNT).min(required_start - 1) {
        for rl in 1..=MAX_REPEAT_LENGTH {
            if required_start - rl * OUTER_MIN_REPEAT_COUNT < 0 {
                continue;
            }
            let Some(repeat) = find_multi_base_repeat(bases, index, rl, MIN_REPEAT_COUNT) else { continue };
            let repeat = extend_repeat_lower(repeat, bases);
            if repeat.count < OUTER_MIN_REPEAT_COUNT {
                continue;
            }
            if !position_within(required_start, repeat.index, repeat.end_index()) {
                continue;
            }
            if max_repeat.as_ref().map_or(true, |m| repeat.total_length() > m.total_length()) {
                max_repeat = Some(repeat);
            }
        }
        index += 1;
    }
    max_repeat
}

fn find_extra_upper_repeat(bases: &[u8], required_index: i32) -> Option<RepeatInfo> {
    let mut max_repeat: Option<RepeatInfo> = None;
    for rl in 1..=MAX_REPEAT_LENGTH {
        let Some(repeat) = find_multi_base_repeat(bases, required_index, rl, OUTER_MIN_REPEAT_COUNT) else {
            continue;
        };
        let repeat = extend_repeat_lower(repeat, bases);
        if max_repeat.as_ref().map_or(true, |m| repeat.total_length() > m.total_length()) {
            max_repeat = Some(repeat);
        }
    }
    max_repeat
}

pub fn find_repeat_boundaries(
    bases: &[u8], required_start: i32, required_end: i32, max_length: i32, min_count: i32,
) -> RepeatBoundaries {
    let mut index = 0.max(required_start - REPEAT_SEARCH_LENGTH);
    let mut all: Option<Vec<RepeatInfo>> = None;
    while index <= (bases.len() as i32 - min_count).min(required_end) {
        for rl in 1..=max_length {
            let Some(repeat) = find_multi_base_repeat(bases, index, rl, min_count) else { continue };
            let repeat = extend_repeat_lower(repeat, bases);
            if !positions_overlap(required_start, required_end, repeat.index, repeat.end_index()) {
                continue;
            }
            match all.as_mut() {
                None => all = Some(vec![repeat]),
                Some(list) => add_if_unique(list, repeat),
            }
        }
        index += 1;
    }

    let mut max_repeat: Option<RepeatInfo> = None;
    let mut lower_index = required_start;
    let mut upper_index = required_end;
    if let Some(list) = all.as_ref() {
        let mut lower: Option<&RepeatInfo> = None;
        let mut upper: Option<&RepeatInfo> = None;
        let mut max: Option<&RepeatInfo> = None;
        for repeat in list {
            if max.map_or(true, |m| {
                repeat.count > m.count || (repeat.count == m.count && repeat.bases.len() > m.bases.len())
            }) {
                max = Some(repeat);
            }
            if position_within(required_start, repeat.index, repeat.end_index())
                && lower.map_or(true, |l| repeat.index < l.index)
            {
                lower = Some(repeat);
            }
            if position_within(required_end, repeat.index, repeat.end_index())
                && upper.map_or(true, |u| repeat.end_index() > u.end_index())
            {
                upper = Some(repeat);
            }
        }
        let lower = lower.or(max).unwrap();
        let upper = upper.or(max).unwrap();
        lower_index = find_post_repeat_index(lower, bases, true);
        upper_index = find_post_repeat_index(upper, bases, false);
        max_repeat = max.cloned();
    }

    let outer_lower_required = if lower_index < required_start { lower_index } else { required_start - 1 };
    if let Some(outer) = find_extra_lower_repeat(bases, outer_lower_required) {
        if outer.index < lower_index {
            lower_index = find_post_repeat_index(&outer, bases, true);
            all.get_or_insert_with(Vec::new).push(outer);
        }
    }
    let outer_upper_required = if upper_index > required_end { upper_index } else { required_end + 1 };
    if let Some(outer) = find_extra_upper_repeat(bases, outer_upper_required) {
        if outer.end_index() > upper_index {
            upper_index = find_post_repeat_index(&outer, bases, false);
            all.get_or_insert_with(Vec::new).push(outer);
        }
    }
    RepeatBoundaries { lower_index, upper_index, max_repeat, all_repeats: all }
}
