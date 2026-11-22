//! SequenceDiffInfo, SequenceDiffType and BaseQualType: how a read differs from the consensus at one point.

use crate::assembly::read_parse_state::ReadParseState;
use crate::assembly::types::repeat_info::RepeatInfo;
use crate::common::fmt::java_fixed;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SequenceDiffType {
    Unset,
    Match,
    Base,
    Delete,
    Insert,
    Repeat,
    Other,
}

impl SequenceDiffType {
    pub fn name(self) -> &'static str {
        match self {
            SequenceDiffType::Unset => "UNSET",
            SequenceDiffType::Match => "MATCH",
            SequenceDiffType::Base => "BASE",
            SequenceDiffType::Delete => "DELETE",
            SequenceDiffType::Insert => "INSERT",
            SequenceDiffType::Repeat => "REPEAT",
            SequenceDiffType::Other => "OTHER",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum BaseQualType {
    High,
    Medium,
    Low,
}

impl BaseQualType {
    pub fn select_higher(first: BaseQualType, second: BaseQualType) -> BaseQualType {
        if first <= second { first } else { second }
    }
    pub fn select_lower(first: BaseQualType, second: BaseQualType) -> BaseQualType {
        if first >= second { first } else { second }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SequenceDiffInfo {
    pub read_index: i32,
    pub consensus_index: i32,
    /// null in Java for a match
    pub bases: Option<String>,
    pub type_: SequenceDiffType,
    pub indel_length: i32,
    pub repeat_count: i32,
    pub repeat_index_begin: i32,
    pub qual_type: BaseQualType,
    pub mismatch_penalty: f64,
}

impl SequenceDiffInfo {
    /// the shared UNSET instance of the Java class; `is_unset` stands for the identity test against it
    pub fn unset() -> SequenceDiffInfo {
        SequenceDiffInfo::new_full(-1, -1, Some(String::new()), SequenceDiffType::Unset, BaseQualType::Low, 0, -1, 0)
    }
    pub fn is_unset(&self) -> bool {
        self.type_ == SequenceDiffType::Unset
    }
    pub fn new(read_index: i32, consensus_index: i32, bases: Option<String>, type_: SequenceDiffType, qual_type: BaseQualType) -> SequenceDiffInfo {
        Self::new_full(read_index, consensus_index, bases, type_, qual_type, 0, -1, 0)
    }
    #[allow(clippy::too_many_arguments)]
    pub fn new_full(
        read_index: i32, consensus_index: i32, bases: Option<String>, type_: SequenceDiffType, qual_type: BaseQualType, repeat_count: i32,
        repeat_index_begin: i32, indel_length: i32,
    ) -> SequenceDiffInfo {
        SequenceDiffInfo { read_index, consensus_index, bases, type_, indel_length, repeat_count, repeat_index_begin, qual_type, mismatch_penalty: 0.0 }
    }
    pub fn from_snv(read: &ReadParseState, consensus_index: i32) -> SequenceDiffInfo {
        Self::new(read.read_index(), consensus_index, Some((read.current_base() as char).to_string()), SequenceDiffType::Base, read.qual_type(read.read_index()))
    }
    pub fn from_match(read_index: i32, consensus_index: i32) -> SequenceDiffInfo {
        Self::new(read_index, consensus_index, None, SequenceDiffType::Match, BaseQualType::High)
    }
    pub fn from_indel(read_index: i32, consensus_index: i32, bases: &str, type_: SequenceDiffType, qual_type: BaseQualType, indel_length: i32) -> SequenceDiffInfo {
        Self::new_full(read_index, consensus_index, Some(bases.to_string()), type_, qual_type, 0, 0, indel_length)
    }
    pub fn bases_str(&self) -> &str {
        self.bases.as_deref().unwrap_or("null")
    }

    pub fn repeat_index(&self, consensus_repeat: &RepeatInfo, build_forwards: bool, get_start: bool) -> i32 {
        if build_forwards == get_start {
            return self.repeat_index_begin;
        }
        if build_forwards {
            self.repeat_index_begin + self.repeat_count * consensus_repeat.repeat_length() - 1
        } else {
            self.repeat_index_begin - self.repeat_count * consensus_repeat.repeat_length() + 1
        }
    }

    pub fn repeat_index_static(read_index: i32, read_repeat_count: i32, consensus_repeat: &RepeatInfo, previous_repeat_bases: i32, build_forwards: bool, get_start: bool) -> i32 {
        if build_forwards == get_start {
            return if build_forwards { read_index - previous_repeat_bases } else { read_index + previous_repeat_bases };
        }
        if build_forwards {
            let repeat_index_start = read_index - previous_repeat_bases;
            repeat_index_start + read_repeat_count * consensus_repeat.repeat_length() - 1
        } else {
            let repeat_index_end = read_index + previous_repeat_bases;
            repeat_index_end - read_repeat_count * consensus_repeat.repeat_length() + 1
        }
    }
}

impl std::fmt::Display for SequenceDiffInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let t = self.type_.name();
        let pen = java_fixed(self.mismatch_penalty, 2);
        match self.type_ {
            SequenceDiffType::Match => write!(f, "index(r={} c={}): {}", self.read_index, self.consensus_index, t),
            SequenceDiffType::Base | SequenceDiffType::Other => write!(f, "index(r={} c={}): {} {} penalty({})", self.read_index, self.consensus_index, t, self.bases_str(), pen),
            SequenceDiffType::Insert | SequenceDiffType::Delete => {
                write!(f, "index(r={} c={}): {}-{} {} penalty({})", self.read_index, self.consensus_index, t, self.indel_length, self.bases_str(), pen)
            }
            SequenceDiffType::Repeat => write!(f, "index(r={} c={}): {} x{} penalty({})", self.read_index, self.consensus_index, t, self.repeat_count, pen),
            SequenceDiffType::Unset => write!(f, "UNSET"),
        }
    }
}
