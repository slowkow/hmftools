//! Process-wide settings that the Java code keeps in static fields (SvConstants.SEQUENCING_TYPE, the static
//! members of AssemblyConfig and the non-final members of AssemblyConstants).

use std::sync::atomic::{AtomicBool, AtomicI32, AtomicU64, AtomicU8, Ordering::Relaxed};
use std::sync::RwLock;
#[cfg(test)]
use std::sync::Mutex;

use crate::common::compat;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SequencingType {
    Illumina,
    Ultima,
    Sbx,
}

impl SequencingType {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "ILLUMINA" => Some(SequencingType::Illumina),
            "ULTIMA" => Some(SequencingType::Ultima),
            "SBX" => Some(SequencingType::Sbx),
            _ => None,
        }
    }
}

static SEQUENCING_TYPE: AtomicU8 = AtomicU8::new(0);

pub fn sequencing_type() -> SequencingType {
    match SEQUENCING_TYPE.load(Relaxed) {
        0 => SequencingType::Illumina,
        1 => SequencingType::Ultima,
        _ => SequencingType::Sbx,
    }
}
pub fn set_sequencing_type(t: SequencingType) {
    SEQUENCING_TYPE.store(
        match t {
            SequencingType::Illumina => 0,
            SequencingType::Ultima => 1,
            SequencingType::Sbx => 2,
        },
        Relaxed,
    );
}
#[inline]
pub fn is_illumina() -> bool {
    SEQUENCING_TYPE.load(Relaxed) == 0
}
#[inline]
pub fn is_ultima() -> bool {
    SEQUENCING_TYPE.load(Relaxed) == 1
}
#[inline]
pub fn is_sbx() -> bool {
    SEQUENCING_TYPE.load(Relaxed) == 2
}
#[inline]
pub fn has_paired_reads() -> bool {
    is_illumina()
}

macro_rules! global_i32 {
    ($name:ident, $get:ident, $set:ident, $init:expr) => {
        static $name: AtomicI32 = AtomicI32::new($init);
        #[inline]
        pub fn $get() -> i32 {
            $name.load(Relaxed)
        }
        pub fn $set(v: i32) {
            $name.store(v, Relaxed)
        }
    };
}
macro_rules! global_bool {
    ($name:ident, $get:ident, $set:ident, $init:expr) => {
        static $name: AtomicBool = AtomicBool::new($init);
        #[inline]
        pub fn $get() -> bool {
            $name.load(Relaxed)
        }
        pub fn $set(v: bool) {
            $name.store(v, Relaxed)
        }
    };
}
macro_rules! global_f64 {
    ($name:ident, $get:ident, $set:ident) => {
        static $name: AtomicU64 = AtomicU64::new(0);
        #[inline]
        pub fn $get() -> f64 {
            f64::from_bits($name.load(Relaxed))
        }
        pub fn $set(v: f64) {
            $name.store(v.to_bits(), Relaxed)
        }
    };
}

// AssemblyConstants (non-final)
global_i32!(MAX_OBSERVED_CONCORDANT_FRAG_LENGTH, max_observed_concordant_frag_length, set_max_observed_concordant_frag_length, 1000);
global_i32!(READ_MISMATCH_MEDIUM_REPEAT_COUNT, read_mismatch_medium_repeat_count, set_read_mismatch_medium_repeat_count, 11);
global_i32!(READ_MISMATCH_LONG_REPEAT_COUNT, read_mismatch_long_repeat_count, set_read_mismatch_long_repeat_count, 16);
global_i32!(READ_MISMATCH_MAX_REPEAT_DIFF, read_mismatch_max_repeat_diff, set_read_mismatch_max_repeat_diff, 5);
global_i32!(PROXIMATE_JUNCTION_DISTANCE, proximate_junction_distance, set_proximate_junction_distance, 50);

// AssemblyApplication.SampleDiscordantRate
global_f64!(SAMPLE_DISCORDANT_RATE, sample_discordant_rate, set_sample_discordant_rate);

// AssemblyConfig statics
global_i32!(PHASE_PROCESSING_LIMIT, phase_processing_limit, set_phase_processing_limit, 0);
global_bool!(WRITE_CANDIDATE_READS, write_candidate_reads, set_write_candidate_reads, false);
global_bool!(ASSEMBLY_BUILD_DEBUG, assembly_build_debug, set_assembly_build_debug, false);
global_bool!(ASSEMBLY_LOG_ALL_READS, assembly_log_all_reads, set_assembly_log_all_reads, false);
global_bool!(DEV_DEBUG, dev_debug, set_dev_debug, false);
global_bool!(PERF_DEBUG, perf_debug, set_perf_debug, false);
global_f64!(PERF_LOG_TIME, perf_log_time, set_perf_log_time);
global_bool!(CLEAR_CACHED_READ, clear_cached_read, set_clear_cached_read, true);

/// SvConstants.maxConcordantFragmentLength
pub fn max_concordant_fragment_length(observed_max_fragment_length: i32) -> i32 {
    observed_max_fragment_length.max(super::constants::MIN_UPPER_FRAGMENT_LENGTH)
}

/// ReadIdTrimmer: removes the leading part of read names that all reads share (up to the third ':').
///
/// Java (finding A5): the prefix comes from the first read name that any thread asks for, and trimming
/// switches itself off at the first name whose third ':' is elsewhere; a name with its third ':' in the
/// same place is cut there whatever comes before it. When tumor and normal, or the lanes of one sample,
/// come from different sequencing runs, ids then depend on processing order (and so on the thread count),
/// and reads of different runs can be given the same id.
///
/// Corrected: the prefix is set once before any read is processed (`set_prefix_from`, from the first
/// record of the tumor BAM), only names that start with that prefix are trimmed, every other name is kept
/// whole, and trimming never switches off. Where all reads share one prefix the ids are Java's.
pub struct ReadIdTrimmer {
    enabled: AtomicBool,
    read_id_trim_index: AtomicI32,
    trimmed_bases: RwLock<String>,
    /// corrected: the number of ':' in a trimmed id (that of the prefix-setting name less three)
    trimmed_delims: AtomicI32,
}

const READ_ID_DELIM: u8 = b':';
const READ_ID_DELIM_COUNT: i32 = 3;

impl ReadIdTrimmer {
    pub const fn new(enabled: bool) -> Self {
        ReadIdTrimmer {
            enabled: AtomicBool::new(enabled),
            read_id_trim_index: AtomicI32::new(-1),
            trimmed_bases: RwLock::new(String::new()),
            trimmed_delims: AtomicI32::new(-1),
        }
    }
    pub fn reset(&self, enabled: bool) {
        self.enabled.store(enabled, Relaxed);
        self.read_id_trim_index.store(-1, Relaxed);
        self.trimmed_bases.write().unwrap().clear();
        self.trimmed_delims.store(-1, Relaxed);
    }
    pub fn enabled(&self) -> bool {
        self.enabled.load(Relaxed)
    }
    /// Corrected (finding A5): fixes the prefix from a read name before any read is processed; a name
    /// without three ':' turns trimming off, as the first name does in Java.
    pub fn set_prefix_from(&self, read_name: &str) {
        if !self.enabled() {
            return;
        }
        let index = find_read_id_trim_index(read_name);
        if index < 0 || index as usize >= read_name.len() {
            self.enabled.store(false, Relaxed);
            return;
        }
        *self.trimmed_bases.write().unwrap() = read_name[..index as usize].to_string();
        let delims = read_name.bytes().filter(|&b| b == READ_ID_DELIM).count() as i32 - READ_ID_DELIM_COUNT;
        self.trimmed_delims.store(delims, Relaxed);
        self.read_id_trim_index.store(index, Relaxed);
    }
    pub fn trim(&self, read_id: &str) -> String {
        if !self.enabled() {
            return read_id.to_string();
        }
        if !compat::java("A5") {
            let prefix = self.trimmed_bases.read().unwrap();
            return match read_id.strip_prefix(prefix.as_str()) {
                Some(rest) if !prefix.is_empty() && !rest.is_empty() => rest.to_string(),
                _ => read_id.to_string(),
            };
        }
        let mut index = self.read_id_trim_index.load(Relaxed);
        if index < 0 {
            let mut prefix = self.trimmed_bases.write().unwrap();
            index = self.read_id_trim_index.load(Relaxed);
            if index < 0 {
                index = find_read_id_trim_index(read_id);
                if index < 0 {
                    self.enabled.store(false, Relaxed);
                    return read_id.to_string();
                }
                *prefix = read_id[..index as usize].to_string();
                self.read_id_trim_index.store(index, Relaxed);
            }
        }
        let i = index as usize;
        if i >= read_id.len() || read_id.as_bytes()[i - 1] != READ_ID_DELIM {
            self.enabled.store(false, Relaxed);
            return read_id.to_string();
        }
        read_id[i..].to_string()
    }
    pub fn restore(&self, trimmed_read_id: &str) -> String {
        if !self.enabled() {
            return trimmed_read_id.to_string();
        }
        let prefix = self.trimmed_bases.read().unwrap();
        if !compat::java("A5") {
            // only a trimmed id gets the prefix back: names of other runs were kept whole (exact where the
            // names have one format, as Illumina names do)
            let delims = trimmed_read_id.bytes().filter(|&b| b == READ_ID_DELIM).count() as i32;
            if prefix.is_empty() || delims != self.trimmed_delims.load(Relaxed) {
                return trimmed_read_id.to_string();
            }
        }
        format!("{}{}", prefix, trimmed_read_id)
    }
}

#[cfg(test)]
mod read_id_trimmer_tests {
    use super::ReadIdTrimmer;
    use crate::common::compat::with_java;

    const RUN_A: &str = "A00260:251:HLYGFDSXY:1:1673:32280:4946";
    const RUN_A_MATE_LANE: &str = "A00260:251:HLYGFDSXY:2:1101:1000:2000";
    // another run, prefix of the same length: Java cuts it at the same place
    const RUN_B: &str = "A00260:252:HLYGFDSXZ:1:1673:32280:4946";
    // another run, prefix one character shorter: Java switches trimming off here
    const RUN_C: &str = "A0026:251:HLYGFDSXY:1:1673:32280:4946";

    #[test]
    fn corrected_trims_only_the_prefix_of_the_first_record() {
        with_java(false, || {
            let t = ReadIdTrimmer::new(true);
            t.set_prefix_from(RUN_A);
            assert_eq!(t.trim(RUN_A), "1:1673:32280:4946");
            assert_eq!(t.trim(RUN_A_MATE_LANE), "2:1101:1000:2000");
            assert_eq!(t.trim(RUN_B), RUN_B);
            assert_eq!(t.trim(RUN_C), RUN_C);
            // still on after names of other runs, and in any order
            assert!(t.enabled());
            assert_eq!(t.trim(RUN_A), "1:1673:32280:4946");
            assert_ne!(t.trim(RUN_A), t.trim(RUN_B));
            assert_eq!(t.restore("1:1673:32280:4946"), RUN_A);
            assert_eq!(t.restore(RUN_B), RUN_B);
            assert_eq!(t.restore(RUN_C), RUN_C);
            assert_eq!(t.restore("null"), "null");
        });
    }

    #[test]
    fn corrected_without_a_prefix_trims_nothing() {
        with_java(false, || {
            let t = ReadIdTrimmer::new(true);
            assert_eq!(t.trim(RUN_A), RUN_A);
            assert_eq!(t.restore(RUN_A), RUN_A);
            let t = ReadIdTrimmer::new(true);
            t.set_prefix_from("SRR123.1");
            assert!(!t.enabled());
            assert_eq!(t.trim(RUN_A), RUN_A);
        });
    }

    #[test]
    fn java_ids_depend_on_order() {
        with_java(true, || {
            // the run of the first name seen decides; a name of the same length from another run is cut
            let t = ReadIdTrimmer::new(true);
            assert_eq!(t.trim(RUN_A), "1:1673:32280:4946");
            assert_eq!(t.trim(RUN_B), "1:1673:32280:4946");
            // a name with its third ':' elsewhere switches trimming off for every later name
            assert_eq!(t.trim(RUN_C), RUN_C);
            assert!(!t.enabled());
            assert_eq!(t.trim(RUN_A), RUN_A);
            // the same names in another order give other ids
            let t = ReadIdTrimmer::new(true);
            assert_eq!(t.trim(RUN_C), "1:1673:32280:4946");
            assert_eq!(t.trim(RUN_A), RUN_A);
        });
    }
}

fn find_read_id_trim_index(read_id: &str) -> i32 {
    let mut delim_count = 0;
    for (i, b) in read_id.bytes().enumerate() {
        if b == READ_ID_DELIM {
            delim_count += 1;
        }
        if delim_count == READ_ID_DELIM_COUNT {
            return i as i32 + 1;
        }
    }
    -1
}

pub static READ_ID_TRIMMER: ReadIdTrimmer = ReadIdTrimmer::new(false);

/// Tests that change the process-wide settings hold this lock.
#[cfg(test)]
pub static TEST_GLOBALS_LOCK: Mutex<()> = Mutex::new(());
