//! PerformanceCounter as used by the assembly tasks (ThreadTask): interval counts and times for -perf_debug,
//! and the "exceeds perf log time" lines of -perf_log_time.

use std::time::Instant;

use crate::assembly::globals::perf_log_time;
use crate::common::fmt::java_fixed;

#[derive(Clone, Debug)]
pub struct PerfCounter {
    name: &'static str,
    interval_count: u64,
    total_time: f64,
    max_time: f64,
}

impl PerfCounter {
    pub fn new(name: &'static str) -> PerfCounter {
        PerfCounter { name, interval_count: 0, total_time: 0.0, max_time: 0.0 }
    }

    /// ThreadTask.stopCheckLog: records the interval started at `start` and logs it if above the perf log time.
    pub fn stop_check_log(&mut self, start: Instant, event_details: impl FnOnce() -> String) {
        let seconds = start.elapsed().as_secs_f64();
        self.max_time = self.max_time.max(seconds);
        self.total_time += seconds;
        self.interval_count += 1;
        let threshold = perf_log_time();
        if threshold > 0.0 && seconds > threshold {
            crate::info!("{}: exceeds perf log time({})", event_details(), java_fixed(seconds, 3));
        }
    }

    pub fn merge(&mut self, other: &PerfCounter) {
        self.interval_count += other.interval_count;
        self.total_time += other.total_time;
        self.max_time = self.max_time.max(other.max_time);
    }

    pub fn log_stats(&self) {
        if self.interval_count == 0 {
            return;
        }
        let avg = self.total_time / self.interval_count as f64;
        crate::info!(
            "PerfStat({}) intervals({}) total({} avg={} max={})",
            self.name,
            self.interval_count,
            java_fixed(self.total_time, 3),
            java_fixed(avg, 3),
            java_fixed(self.max_time, 3)
        );
    }
}
