//! ReadStats: counts of reads and read adjustments made while loading.

#[derive(Default, Clone, Debug)]
pub struct ReadStats {
    pub total_reads: i32,
    pub low_base_qual_filtered: i32,
    pub poly_g_trimmed: i32,
    pub low_base_qual_trimmed: i32,
    /// reads with adapter bases trimmed (only counted when finding A12 is corrected)
    pub adapter_trimmed: i32,
    pub indel_soft_clip_converted: i32,
    pub decoy_sequences: i32,
    pub ref_base_alignment_fails: i32,
    pub identical_supplementaries: i32,
}

impl ReadStats {
    pub fn merge(&mut self, other: &ReadStats) {
        self.total_reads += other.total_reads;
        self.low_base_qual_filtered += other.low_base_qual_filtered;
        self.poly_g_trimmed += other.poly_g_trimmed;
        self.low_base_qual_trimmed += other.low_base_qual_trimmed;
        self.adapter_trimmed += other.adapter_trimmed;
        self.indel_soft_clip_converted += other.indel_soft_clip_converted;
        self.decoy_sequences += other.decoy_sequences;
        self.ref_base_alignment_fails += other.ref_base_alignment_fails;
        self.identical_supplementaries += other.identical_supplementaries;
    }
}

impl std::fmt::Display for ReadStats {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // Java has no adapter count: its adapter trim never happens (finding A12)
        let adapter = if crate::common::compat::java("A12") { String::new() } else { format!(" adapter={}", self.adapter_trimmed) };
        write!(
            f,
            "reads({}) lowBaseQual({}) trim(polyG={} lowBase={}{adapter}) indelSoftClip({}) decoySequences({}) refBaseAlignFails({})",
            self.total_reads,
            self.low_base_qual_filtered,
            self.poly_g_trimmed,
            self.low_base_qual_trimmed,
            self.indel_soft_clip_converted,
            self.decoy_sequences,
            self.ref_base_alignment_fails
        )
    }
}
