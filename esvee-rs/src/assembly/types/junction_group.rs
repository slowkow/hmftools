//! JunctionGroup: junctions close enough to share one BAM slice.

use crate::assembly::common_utils::positions_overlap;
use crate::assembly::constants::BAM_READ_JUNCTION_BUFFER;
use crate::assembly::ptr::P;
use crate::assembly::read::read::ReadP;
use crate::assembly::types::junction::JuncP;
use crate::assembly::types::junction_assembly::AsmP;
use crate::assembly::types::remote_region::RemoteRegion;

pub type JuncGroupP = P<JunctionGroup>;

pub struct JunctionGroup {
    start: i32,
    end: i32,
    junctions: Vec<JuncP>,
    candidate_reads: Vec<ReadP>,
    final_candidate_read_count: usize,
    index: i32,
    junction_assemblies: Vec<AsmP>,
}

impl JunctionGroup {
    pub fn new(junction: JuncP) -> JunctionGroup {
        JunctionGroup {
            start: junction.position,
            end: junction.position,
            junctions: vec![junction],
            candidate_reads: Vec::new(),
            final_candidate_read_count: 0,
            index: -1,
            junction_assemblies: Vec::new(),
        }
    }
    pub fn start(&self) -> i32 {
        self.start
    }
    pub fn end(&self) -> i32 {
        self.end
    }
    pub fn junctions(&self) -> &Vec<JuncP> {
        &self.junctions
    }
    pub fn count(&self) -> usize {
        self.junctions.len()
    }
    pub fn chromosome(&self) -> &'static str {
        self.junctions[0].chromosome
    }
    pub fn min_position(&self) -> i32 {
        self.start
    }
    pub fn max_position(&self) -> i32 {
        self.end
    }
    pub fn range(&self) -> i32 {
        self.end - self.start + 1
    }
    pub fn read_range_start(&self) -> i32 {
        self.start - BAM_READ_JUNCTION_BUFFER
    }
    pub fn read_range_end(&self) -> i32 {
        self.end + BAM_READ_JUNCTION_BUFFER
    }
    pub fn set_index(&mut self, index: i32) {
        self.index = index;
    }
    pub fn index(&self) -> i32 {
        self.index
    }
    pub fn add_junction(&mut self, junction: JuncP) {
        self.end = self.end.max(junction.position);
        self.start = self.start.min(junction.position);
        self.junctions.push(junction);
    }
    pub fn candidate_reads(&self) -> &Vec<ReadP> {
        &self.candidate_reads
    }
    pub fn add_candidate_read(&mut self, read: ReadP) {
        self.candidate_reads.push(read);
    }
    pub fn clear_candidate_reads(&mut self) {
        self.final_candidate_read_count = self.candidate_reads.len();
        self.candidate_reads = Vec::new();
    }
    pub fn candidate_read_count(&self) -> usize {
        self.final_candidate_read_count.max(self.candidate_reads.len())
    }
    pub fn add_junction_assemblies(&mut self, assemblies: &[AsmP]) {
        self.junction_assemblies.extend_from_slice(assemblies);
    }
    pub fn junction_assemblies(&self) -> &Vec<AsmP> {
        &self.junction_assemblies
    }
    pub fn overlaps_remote_region(&self, region: &RemoteRegion) -> bool {
        positions_overlap(self.read_range_start(), self.read_range_end(), region.start(), region.end())
    }

    /// Index of the last group with a start at or before `read_start`, or zero if all start later.
    /// (Collections.binarySearch on the start positions; with equal starts the entry found is the one the
    /// Java search lands on.)
    pub fn binary_search(read_start: i32, regions: &[JuncGroupP]) -> usize {
        let mut low = 0isize;
        let mut high = regions.len() as isize - 1;
        while low <= high {
            let mid = ((low + high) as usize >> 1) as isize;
            let mid_val = regions[mid as usize].start();
            if mid_val < read_start {
                low = mid + 1;
            } else if mid_val > read_start {
                high = mid - 1;
            } else {
                return mid as usize;
            }
        }
        let insertion_index = low;
        if insertion_index == 0 {
            return 0;
        }
        (insertion_index - 1) as usize
    }

    pub fn build_junction_groups(junctions: &[JuncP], max_distance: i32) -> Vec<JuncGroupP> {
        let mut groups: Vec<JuncGroupP> = Vec::new();
        let mut current_group = P::new(JunctionGroup::new(junctions[0]));
        groups.push(current_group);
        for junction in junctions.iter().skip(1) {
            if junction.position - current_group.max_position() <= max_distance {
                current_group.m().add_junction(*junction);
            } else {
                current_group = P::new(JunctionGroup::new(*junction));
                groups.push(current_group);
            }
        }
        groups
    }
}

impl std::fmt::Display for JunctionGroup {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}:{}-{} range({}) count({}) assemblies({}) reads({})",
            self.chromosome(),
            self.start,
            self.end,
            self.range(),
            self.junctions.len(),
            self.junction_assemblies.len(),
            self.candidate_read_count()
        )
    }
}
