//! Log lines in the layout of the Java tools: `HH:MM:SS.mmm [INFO ] message`, on stdout.

use std::sync::atomic::{AtomicBool, Ordering};

pub static DEBUG: AtomicBool = AtomicBool::new(false);

pub fn debug_enabled() -> bool {
    DEBUG.load(Ordering::Relaxed)
}

pub fn write(level: &str, msg: std::fmt::Arguments) {
    use std::io::Write;
    let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap();
    let secs = now.as_secs();
    let mut tm: libc_tm::Tm = libc_tm::local(secs as i64);
    tm.millis = now.subsec_millis();
    let out = std::io::stdout();
    let mut out = out.lock();
    let _ = writeln!(out, "{:02}:{:02}:{:02}.{:03} [{}] {}", tm.hour, tm.min, tm.sec, tm.millis, level, msg);
}

mod libc_tm {
    pub struct Tm {
        pub hour: i32,
        pub min: i32,
        pub sec: i32,
        pub millis: u32,
    }
    unsafe extern "C" {
        fn localtime_r(t: *const i64, tm: *mut [i32; 16]) -> *mut [i32; 16];
    }
    pub fn local(secs: i64) -> Tm {
        let mut raw = [0i32; 16];
        unsafe { localtime_r(&secs, &mut raw) };
        Tm { sec: raw[0], min: raw[1], hour: raw[2], millis: 0 }
    }
}

#[macro_export]
macro_rules! info {
    ($($arg:tt)*) => { $crate::log::write("INFO ", format_args!($($arg)*)) };
}
#[macro_export]
macro_rules! warn {
    ($($arg:tt)*) => { $crate::log::write("WARN ", format_args!($($arg)*)) };
}
#[macro_export]
macro_rules! error {
    ($($arg:tt)*) => { $crate::log::write("ERROR", format_args!($($arg)*)) };
}
#[macro_export]
macro_rules! debug {
    ($($arg:tt)*) => { if $crate::log::debug_enabled() { $crate::log::write("DEBUG", format_args!($($arg)*)) } };
}
