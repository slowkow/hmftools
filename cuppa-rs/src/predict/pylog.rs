//! Log lines as pycuppa's logging setup writes them (cuppa/logger.py): to stdout, and to a file with
//! `--log_to_file` / `--log_path`, in the `--log_format` layout (default
//! `%(asctime)s [%(threadName)s] %(levelname)s %(name)s | %(message)s`, times as %H:%M:%S).

use std::io::Write;
use std::sync::Mutex;

pub const DEFAULT_FORMAT: &str = "%(asctime)s [%(threadName)s] %(levelname)s %(name)s | %(message)s";

struct Setup {
    format: String,
    file: Option<std::fs::File>,
}

static SETUP: Mutex<Option<Setup>> = Mutex::new(None);

/// initialize_logging(filename, format): the file is truncated (mode 'w')
pub fn init(format: &str, file: Option<&str>) -> std::io::Result<()> {
    let f = match file {
        Some(p) => Some(std::fs::File::create(p)?),
        None => None,
    };
    *SETUP.lock().unwrap() = Some(Setup { format: format.to_string(), file: f });
    Ok(())
}

fn now_parts() -> (i32, i32, i32, u32) {
    let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap();
    unsafe extern "C" {
        fn localtime_r(t: *const i64, tm: *mut [i32; 16]) -> *mut [i32; 16];
    }
    let secs = now.as_secs() as i64;
    let mut raw = [0i32; 16];
    unsafe { localtime_r(&secs, &mut raw) };
    (raw[2], raw[1], raw[0], now.subsec_millis())
}

/// %-style substitution of the logging record's attributes
fn render(format: &str, level: &str, name: &str, func: &str, msg: &str) -> String {
    let (h, m, s, ms) = now_parts();
    let mut out = String::new();
    let b = format.as_bytes();
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'%' && i + 1 < b.len() && b[i + 1] == b'%' {
            out.push('%');
            i += 2;
            continue;
        }
        if b[i] == b'%' && i + 1 < b.len() && b[i + 1] == b'(' {
            if let Some(close) = format[i + 2..].find(')') {
                let key = &format[i + 2..i + 2 + close];
                let mut j = i + 2 + close + 1;
                let spec_start = j;
                while j < b.len() && !b[j].is_ascii_alphabetic() {
                    j += 1;
                }
                if j < b.len() {
                    let spec = &format[spec_start..j];
                    let conv = b[j] as char;
                    let value = match key {
                        "asctime" => format!("{h:02}:{m:02}:{s:02}"),
                        "msecs" => ms.to_string(),
                        "threadName" => "MainThread".to_string(),
                        "levelname" => level.to_string(),
                        "name" => name.to_string(),
                        "funcName" => func.to_string(),
                        "message" => msg.to_string(),
                        "process" => std::process::id().to_string(),
                        _ => String::new(),
                    };
                    let width: usize = spec.trim_start_matches('0').trim_start_matches('-').parse().unwrap_or(0);
                    let v = if conv == 'd' && spec.starts_with('0') && value.len() < width {
                        format!("{}{}", "0".repeat(width - value.len()), value)
                    } else if value.len() < width {
                        if spec.starts_with('-') { format!("{value:<width$}") } else { format!("{value:>width$}") }
                    } else {
                        value
                    };
                    out.push_str(&v);
                    i = j + 1;
                    continue;
                }
            }
        }
        let ch = format[i..].chars().next().unwrap();
        out.push(ch);
        i += ch.len_utf8();
    }
    out
}

fn write(level: &str, name: &str, func: &str, msg: &str) {
    let mut guard = SETUP.lock().unwrap();
    let fmt = guard.as_ref().map(|s| s.format.clone()).unwrap_or_else(|| DEFAULT_FORMAT.to_string());
    let line = render(&fmt, level, name, func, msg);
    {
        let out = std::io::stdout();
        let mut out = out.lock();
        let _ = writeln!(out, "{line}");
    }
    if let Some(Setup { file: Some(f), .. }) = guard.as_mut() {
        let _ = writeln!(f, "{line}");
    }
}

pub fn info(name: &str, func: &str, msg: &str) {
    write("[INFO ]", name, func, msg);
}
pub fn debug(name: &str, func: &str, msg: &str) {
    write("[DEBUG]", name, func, msg);
}
pub fn warning(name: &str, func: &str, msg: &str) {
    write("[WARN ]", name, func, msg);
}
pub fn error(name: &str, func: &str, msg: &str) {
    write("[ERROR]", name, func, msg);
}

#[cfg(test)]
mod tests {
    #[test]
    fn layout() {
        let s = super::render("%(levelname)s %(name)s.%(funcName)s | %(message)s %%", "[INFO ]", "A", "f", "m");
        assert_eq!(s, "[INFO ] A.f | m %");
        let t = super::render("%(msecs)03d", "", "", "", "");
        assert_eq!(t.len(), 3);
    }
}
