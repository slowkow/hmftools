//! JunctionSequence: an assembly's bases, optionally reverse-complemented, with the sub-sequence used to
//! match it against another assembly.

use std::cell::OnceCell;

use crate::assembly::assembly_utils::{reverse_complement_bases, reverse_complement_str};
use crate::assembly::orientation::Orientation;
use crate::assembly::types::junction_assembly::{bytes_to_string, JunctionAssembly};
use crate::assembly::types::repeat_info::{find_repeats, RepeatInfo};

pub const PHASED_ASSEMBLY_MATCH_SEQ_LENGTH: i32 = 100;

#[derive(Clone, Copy, PartialEq, Eq)]
enum MatchSequenceMode {
    Straddle,
    FullExtension,
    OuterExtension,
    InnerExtension,
}

pub struct JunctionSequence {
    pub reversed: bool,
    pub junc_orient: Orientation,
    pub full_sequence: String,
    pub extension_length: i32,
    pub ref_base_length: i32,
    pub base_length: i32,
    original_bases: Vec<u8>,
    original_base_quals: Vec<u8>,
    junction_index: i32,
    match_seq_index_start: i32,
    match_seq_index_end: i32,
    repeat_info: OnceCell<Vec<RepeatInfo>>,
    bases: OnceCell<Vec<u8>>,
    base_quals: OnceCell<Vec<u8>>,
}

impl JunctionSequence {
    pub fn form_outer_extension_match_sequence(assembly: &JunctionAssembly, reverse_compliment: bool) -> JunctionSequence {
        Self::from_assembly(assembly, reverse_compliment, 0, PHASED_ASSEMBLY_MATCH_SEQ_LENGTH, MatchSequenceMode::OuterExtension)
    }
    pub fn form_straddling_match_sequence(assembly: &JunctionAssembly) -> JunctionSequence {
        Self::from_assembly(assembly, false, PHASED_ASSEMBLY_MATCH_SEQ_LENGTH / 2, -1, MatchSequenceMode::Straddle)
    }
    pub fn form_straddling_match_sequence_ext(assembly: &JunctionAssembly, reverse_compliment: bool, max_match_seq_ref_base_length: i32, max_match_seq_extension_length: i32) -> JunctionSequence {
        Self::from_assembly(assembly, reverse_compliment, max_match_seq_ref_base_length, max_match_seq_extension_length, MatchSequenceMode::Straddle)
    }
    pub fn form_full_extension_match_sequence(assembly: &JunctionAssembly, reverse_compliment: bool) -> JunctionSequence {
        Self::from_assembly(assembly, reverse_compliment, 0, -1, MatchSequenceMode::FullExtension)
    }
    pub fn form_inner_extension_match_sequence(assembly: &JunctionAssembly, reverse_compliment: bool) -> JunctionSequence {
        Self::from_assembly(assembly, reverse_compliment, 0, PHASED_ASSEMBLY_MATCH_SEQ_LENGTH, MatchSequenceMode::InnerExtension)
    }

    fn from_assembly(assembly: &JunctionAssembly, reverse_compliment: bool, max_match_seq_ref_base_length: i32, max_match_seq_extension_length: i32, mode: MatchSequenceMode) -> JunctionSequence {
        let reversed = reverse_compliment;
        let ref_base_length = assembly.ref_base_length();
        let extension_length = assembly.extension_length();
        let base_length = ref_base_length + extension_length;
        let junction_index = assembly.junction_index();
        let full_sequence = if !reversed { assembly.form_full_sequence() } else { reverse_complement_str(&assembly.form_full_sequence()) };

        let match_seq_ext_length;
        let mut match_seq_ref_extension = 0;
        if max_match_seq_extension_length > 0 {
            if max_match_seq_extension_length > extension_length {
                match_seq_ext_length = extension_length;
                match_seq_ref_extension = max_match_seq_extension_length - extension_length;
            } else {
                match_seq_ext_length = max_match_seq_extension_length;
            }
        } else {
            match_seq_ext_length = extension_length;
        }
        let match_seq_ref_length = ref_base_length.min(max_match_seq_ref_base_length + match_seq_ref_extension);
        let match_index_start;
        let match_index_end;
        let fwd = assembly.is_forward_junction();
        if mode == MatchSequenceMode::OuterExtension && match_seq_ext_length < extension_length {
            if fwd {
                match_index_end = base_length - 1;
                match_index_start = match_index_end - match_seq_ext_length + 1;
            } else {
                match_index_start = 0;
                match_index_end = match_index_start + match_seq_ext_length - 1;
            }
        } else if mode == MatchSequenceMode::InnerExtension && match_seq_ext_length < extension_length {
            if fwd {
                match_index_start = junction_index + 1;
                match_index_end = match_index_start + match_seq_ext_length - 1;
            } else {
                match_index_end = junction_index - 1;
                match_index_start = match_index_end - match_seq_ext_length + 1;
            }
        } else if mode == MatchSequenceMode::FullExtension {
            if fwd {
                match_index_start = junction_index + 1;
                match_index_end = base_length - 1;
            } else {
                match_index_start = 0;
                match_index_end = junction_index - 1;
            }
        } else if fwd {
            match_index_start = junction_index - match_seq_ref_length + 1;
            match_index_end = junction_index + match_seq_ext_length;
        } else {
            match_index_start = junction_index - match_seq_ext_length;
            match_index_end = junction_index + match_seq_ref_length - 1;
        }
        let mut js = JunctionSequence {
            reversed,
            junc_orient: assembly.junction().orient,
            full_sequence,
            extension_length,
            ref_base_length,
            base_length,
            original_bases: assembly.bases().to_vec(),
            original_base_quals: assembly.base_quals().to_vec(),
            junction_index,
            match_seq_index_start: 0,
            match_seq_index_end: 0,
            repeat_info: OnceCell::new(),
            bases: OnceCell::new(),
            base_quals: OnceCell::new(),
        };
        if !reversed {
            js.match_seq_index_start = match_index_start;
            js.match_seq_index_end = match_index_end;
        } else {
            js.match_seq_index_start = js.index_reversed(match_index_end);
            js.match_seq_index_end = js.index_reversed(match_index_start);
        }
        js
    }

    pub fn from_bases(bases: &[u8], base_quals: &[u8], orientation: Orientation, reverse_compliment: bool) -> JunctionSequence {
        let full_sequence = if !reverse_compliment { bytes_to_string(bases) } else { bytes_to_string(&reverse_complement_bases(bases)) };
        JunctionSequence {
            reversed: reverse_compliment,
            junc_orient: orientation,
            full_sequence,
            extension_length: bases.len() as i32,
            ref_base_length: 0,
            base_length: bases.len() as i32,
            original_bases: bases.to_vec(),
            original_base_quals: base_quals.to_vec(),
            junction_index: if orientation.is_reverse() { 0 } else { bases.len() as i32 - 1 },
            match_seq_index_start: 0,
            match_seq_index_end: 0,
            repeat_info: OnceCell::new(),
            bases: OnceCell::new(),
            base_quals: OnceCell::new(),
        }
    }

    pub fn junction_index(&self) -> i32 {
        if !self.reversed { self.junction_index } else { self.base_length - self.junction_index - 1 }
    }
    pub fn match_sequence(&self) -> &str {
        &self.full_sequence[self.match_seq_index_start as usize..(self.match_seq_index_end + 1) as usize]
    }
    pub fn match_seq_start_index(&self) -> i32 {
        self.match_seq_index_start
    }
    pub fn match_seq_end_index(&self) -> i32 {
        self.match_seq_index_end
    }
    pub fn bases(&self) -> &[u8] {
        if !self.reversed {
            return &self.original_bases;
        }
        self.bases.get_or_init(|| self.full_sequence.as_bytes().to_vec())
    }
    pub fn base_quals(&self) -> &[u8] {
        if !self.reversed {
            return &self.original_base_quals;
        }
        self.base_quals.get_or_init(|| self.original_base_quals.iter().rev().copied().collect())
    }
    pub fn repeat_info(&self) -> &Vec<RepeatInfo> {
        self.repeat_info.get_or_init(|| find_repeats(self.full_sequence.as_bytes()))
    }
    fn index_reversed(&self, index: i32) -> i32 {
        let junc_index_diff = index - self.junction_index;
        self.junction_index() - junc_index_diff
    }
    pub fn index_reverted(&self, index: i32) -> i32 {
        let junc_index_diff = index - self.junction_index();
        self.junction_index - junc_index_diff
    }
}

impl std::fmt::Display for JunctionSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "len({} ref={} ext={} juncIndex={}) {} matchSeq({} - {})",
            self.base_length,
            self.ref_base_length,
            self.extension_length,
            self.junction_index,
            if self.reversed { "rev" } else { "fwd" },
            self.match_seq_start_index(),
            self.match_seq_end_index()
        )
    }
}
