//! DepthTracker: aligned bases per fixed window of a partition.

#[cfg(test)]
use super::rec::*;
#[cfg(test)]
use rust_htslib::bam;

pub const DEPTH_WINDOW_SIZE: i32 = 200;

pub struct DepthTracker {
    start: i32,
    end: i32,
    window_size: i32,
    window_depth: Vec<i64>,
}

impl DepthTracker {
    pub fn new(start: i32, end: i32, window_size: i32) -> Self {
        let length = end - start + 1;
        let windows = if crate::common::compat::java("P3") {
            // Java: (int)ceil(length / size) divides integers first, so a last partial window is left out
            length / window_size
        } else {
            (length + window_size - 1) / window_size
        };
        DepthTracker { start, end, window_size, window_depth: vec![0; windows.max(0) as usize] }
    }

    fn window_index(&self, position: i32) -> Option<usize> {
        if position < self.start || position > self.end {
            return None;
        }
        Some(((position - self.start) / self.window_size) as usize)
    }

    pub fn calc_depth(&self, position: i32) -> f64 {
        match self.window_index(position) {
            Some(i) if i < self.window_depth.len() => {
                // a last partial window (only present when P3 is corrected) is averaged over its own length
                let window_start = self.start + i as i32 * self.window_size;
                let window_length = (window_start + self.window_size - 1).min(self.end) - window_start + 1;
                self.window_depth[i] as f64 / window_length as f64
            }
            _ => 0.0,
        }
    }

    #[cfg(test)]
    pub fn process_read(&mut self, rec: &bam::Record) {
        if has(rec, FLAG_UNMAPPED) {
            return;
        }
        self.process_range(start(rec), end(rec));
    }

    pub fn process_range(&mut self, read_start: i32, read_end: i32) {
        let Some(mut index) = self.window_index(read_start) else { return };
        if index >= self.window_depth.len() {
            return;
        }
        let mut window_start = self.start + index as i32 * self.window_size;
        // the clamp only acts on a last partial window
        let mut window_end = (window_start + self.window_size - 1).min(self.end);
        loop {
            let window_range = window_end.min(read_end) - read_start.max(window_start) + 1;
            if window_range < 0 || window_range > self.window_size {
                crate::error!("depth range failed: read({read_start}-{read_end}) window(idx={index} {window_start}-{window_end}) read range({window_range})");
                return;
            }
            self.window_depth[index] += window_range as i64;
            if read_end <= window_end {
                break;
            }
            index += 1;
            if index >= self.window_depth.len() {
                break;
            }
            window_start = window_end + 1;
            window_end = (window_start + self.window_size - 1).min(self.end);
        }
    }
}
