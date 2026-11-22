//! ReadAssemblyIndices: where a read's bases fall in an assembly's base array.

use crate::assembly::read::read::Read;
use crate::assembly::read::read_utils::INVALID_INDEX;
use crate::assembly::types::junction_assembly::JunctionAssembly;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadAssemblyIndices {
    pub read_index_start: i32,
    pub read_index_end: i32,
    pub assembly_index_start: i32,
    pub junction_index: i32,
}

pub const INVALID_INDICES: ReadAssemblyIndices =
    ReadAssemblyIndices { read_index_start: INVALID_INDEX, read_index_end: INVALID_INDEX, assembly_index_start: INVALID_INDEX, junction_index: INVALID_INDEX };

impl ReadAssemblyIndices {
    pub fn new(read_index_start: i32, read_index_end: i32, assembly_index_start: i32) -> Self {
        ReadAssemblyIndices { read_index_start, read_index_end, assembly_index_start, junction_index: INVALID_INDEX }
    }
    /// as in Java, true when the start index is the invalid value
    pub fn is_valid(&self) -> bool {
        self.read_index_start == INVALID_INDEX
    }
    /// the Java identity test against INVALID_INDICES
    pub fn is_invalid_indices(&self) -> bool {
        *self == INVALID_INDICES
    }
    pub fn junction_read_start_distance(&self, assembly_junction_index: i32) -> i32 {
        assembly_junction_index - self.assembly_index_start + self.read_index_start
    }

    pub fn get_ref_read_indices(assembly: &JunctionAssembly, ref_base_position: i32, read: &Read) -> Option<ReadAssemblyIndices> {
        let junction_position = assembly.junction().position;
        let alignment_start = read.alignment_start();
        let left_soft_clip_offset = read.left_clip_length();
        let right_soft_clip_offset = read.right_clip_length();
        let read_index_end = read.bases_length() - right_soft_clip_offset - 1;
        if assembly.is_forward_junction() {
            if alignment_start >= junction_position {
                return None;
            }
            if alignment_start < ref_base_position {
                let read_index = ref_base_position - alignment_start + left_soft_clip_offset;
                return Some(ReadAssemblyIndices::new(read_index, read_index_end, 0));
            }
            Some(ReadAssemblyIndices::new(left_soft_clip_offset, read_index_end, alignment_start - ref_base_position))
        } else {
            if read.alignment_end() <= junction_position {
                return None;
            }
            if alignment_start < junction_position {
                let read_index = junction_position - alignment_start + left_soft_clip_offset;
                return Some(ReadAssemblyIndices::new(read_index, read_index_end, assembly.junction_index()));
            }
            let assembly_index = assembly.junction_index() + alignment_start - junction_position;
            Some(ReadAssemblyIndices::new(left_soft_clip_offset, read_index_end, assembly_index))
        }
    }
}

impl std::fmt::Display for ReadAssemblyIndices {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "read({}-{}) assembly({} junc={})", self.read_index_start, self.read_index_end, self.assembly_index_start, self.junction_index)
    }
}
