//! hmf-common ConfigBuilder: options registered by the program, parsed from `-name value` and `-flag`
//! as Java parses them (a value may not start with '-' unless the item is a number), with the required
//! items, path checks (wildcards and '$' accepted), integer and decimal checks, `-version` and `-help`
//! as the sole argument, and the log level options (ConfigUtils.addLoggingOptions / setLogLevel).

use anyhow::{anyhow, Result};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ItemType {
    Str,
    Path,
    Paths,
    Integer,
    Decimal,
    Flag,
}

impl ItemType {
    /// ConfigItemType names as ConfigBuilder prints them
    pub fn name(self) -> &'static str {
        match self {
            ItemType::Str => "STRING",
            ItemType::Path => "PATH",
            ItemType::Paths => "PATHS",
            ItemType::Integer => "INTEGER",
            ItemType::Decimal => "DECIMAL",
            ItemType::Flag => "FLAG",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Item {
    pub name: String,
    pub ty: ItemType,
    pub required: bool,
    pub desc: String,
    pub default: Option<String>,
    pub path_prefix: Option<String>,
    value: Option<String>,
}

#[derive(Default)]
pub struct ConfigBuilder {
    app_name: String,
    version: String,
    items: Vec<Item>,
}

pub enum Parsed {
    Run,
    /// -version or -help was printed; exit 0
    Exit,
}

/// Java's String.valueOf(double) for a registered decimal default
fn double_default(d: f64) -> String {
    crate::fmt::java_double_to_string(d)
}

impl ConfigBuilder {
    pub fn new(app_name: &str, version: &str) -> ConfigBuilder {
        ConfigBuilder { app_name: app_name.to_string(), version: version.to_string(), items: Vec::new() }
    }

    /// addConfigItem: a second registration of a name is ignored
    pub fn add(&mut self, ty: ItemType, name: &str, required: bool, desc: &str, default: Option<String>) {
        if self.items.iter().any(|i| i.name == name) {
            return;
        }
        self.items.push(Item { name: name.to_string(), ty, required, desc: desc.to_string(), default, path_prefix: None, value: None });
    }
    pub fn add_str(&mut self, name: &str, required: bool, desc: &str) {
        self.add(ItemType::Str, name, required, desc, None);
    }
    pub fn add_str_default(&mut self, name: &str, required: bool, desc: &str, default: &str) {
        self.add(ItemType::Str, name, required, desc, Some(default.to_string()));
    }
    pub fn add_path(&mut self, name: &str, required: bool, desc: &str) {
        self.add(ItemType::Path, name, required, desc, None);
    }
    pub fn add_paths(&mut self, name: &str, required: bool, desc: &str) {
        self.add(ItemType::Paths, name, required, desc, None);
    }
    pub fn add_prefixed_path(&mut self, name: &str, required: bool, desc: &str, prefix: &str) {
        self.add(ItemType::Path, name, required, desc, None);
        if let Some(i) = self.items.iter_mut().find(|i| i.name == name) {
            i.path_prefix = Some(prefix.to_string());
        }
    }
    pub fn add_flag(&mut self, name: &str, desc: &str) {
        self.add(ItemType::Flag, name, false, desc, Some("false".to_string()));
    }
    pub fn add_integer(&mut self, name: &str, desc: &str, default: i32) {
        self.add(ItemType::Integer, name, false, desc, Some(default.to_string()));
    }
    pub fn add_decimal(&mut self, name: &str, desc: &str, default: f64) {
        self.add(ItemType::Decimal, name, false, desc, Some(double_default(default)));
    }
    pub fn add_required_integer(&mut self, name: &str, desc: &str) {
        self.add(ItemType::Integer, name, true, desc, None);
    }
    pub fn add_required_decimal(&mut self, name: &str, desc: &str) {
        self.add(ItemType::Decimal, name, true, desc, None);
    }

    /// ConfigUtils.addLoggingOptions
    pub fn add_logging_options(&mut self) {
        self.add_flag("log_debug", "Log at DEBUG level");
        self.add_str("log_level", false, "Specify log level: ERROR, WARN, INFO, DEBUG or TRACE");
    }

    /// TaskExecutor.addThreadOptions
    pub fn add_thread_options(&mut self, default: i32) {
        self.add_integer("threads", "Number of threads, 0 or 1 not multi-threaded", default);
    }

    pub fn is_registered(&self, name: &str) -> bool {
        self.items.iter().any(|i| i.name == name)
    }

    fn item(&self, name: &str) -> &Item {
        self.items.iter().find(|i| i.name == name).unwrap_or_else(|| panic!("invalid config item({name}) requested"))
    }

    /// getValue: the value given, else the default
    pub fn value(&self, name: &str) -> Option<String> {
        let item = self.item(name);
        item.value.clone().or_else(|| item.default.clone())
    }
    pub fn value_or(&self, name: &str, default: &str) -> String {
        self.value(name).unwrap_or_else(|| default.to_string())
    }
    /// hasValue: whether the option was on the command line (a default does not count)
    pub fn has_value(&self, name: &str) -> bool {
        self.item(name).value.is_some()
    }
    pub fn has_flag(&self, name: &str) -> bool {
        self.value(name).is_some_and(|v| v.eq_ignore_ascii_case("true"))
    }
    pub fn integer(&self, name: &str) -> i32 {
        self.value(name).and_then(|v| crate::java_text::parse_int(&v)).unwrap_or(0)
    }
    pub fn decimal(&self, name: &str) -> f64 {
        self.value(name).and_then(|v| crate::java_text::parse_double(&v)).unwrap_or(0.0)
    }
    /// ConfigBuilder.getConfigDecimal
    pub fn decimal_or(&self, name: &str, default: f64) -> f64 {
        if self.has_value(name) { self.decimal(name) } else { default }
    }
    /// ConfigBuilder.getConfigInteger
    pub fn integer_or(&self, name: &str, default: i32) -> i32 {
        if self.has_value(name) { self.integer(name) } else { default }
    }
    /// sets a value as the command line would
    pub fn set_value(&mut self, name: &str, value: &str) {
        if let Some(i) = self.items.iter_mut().find(|i| i.name == name) {
            i.value = Some(value.to_string());
        }
    }

    /// the item list as printItems writes it
    pub fn item_lines(&self) -> Vec<String> {
        let mut out: Vec<String> = self
            .items
            .iter()
            .map(|item| {
                let mut s = format!(" -{}:", item.name);
                if item.ty != ItemType::Str {
                    s.push_str(&format!(" type({})", item.ty.name()));
                }
                if let Some(d) = &item.default {
                    if item.ty != ItemType::Flag {
                        s.push_str(&format!(" default({d})"));
                    }
                }
                if item.required {
                    s.push_str(" REQUIRED");
                }
                s.push_str(&format!(" desc({})", item.desc));
                s
            })
            .collect();
        out.sort();
        out.insert(0, "registered config items:".into());
        out
    }

    fn valid_path(&self, item: &Item, value: &str) -> bool {
        let paths: Vec<&str> = if item.ty == ItemType::Paths { value.split(',').collect() } else { vec![value] };
        let prefix = item.path_prefix.as_ref().and_then(|p| self.value(p));
        paths.iter().all(|p| {
            let path = match &prefix {
                Some(pre) => format!("{}{}", check_add_dir_separator(pre), p),
                None => p.to_string(),
            };
            path.contains('*') || path.contains('$') || std::path::Path::new(&path).exists()
        })
    }

    /// checkAndParseCommandLine: Err (already logged) for what Java rejects with exit status 1
    pub fn parse(&mut self, argv: &[String]) -> Result<Parsed> {
        if argv.len() == 1 && (argv[0] == "-version" || argv[0] == "-help") {
            println!("{} version {}", self.app_name, self.version);
            if argv[0] == "-help" {
                for line in self.item_lines() {
                    println!("{line}");
                }
            }
            return Ok(Parsed::Exit);
        }
        let mut matched: Option<usize> = None;
        let mut incorrect = false;
        for arg in argv {
            if let Some(i) = matched {
                let ty = self.items[i].ty;
                if arg.starts_with('-') && ty != ItemType::Decimal && ty != ItemType::Integer {
                    crate::error!("config item({}) has invalid argument: {arg}", self.items[i].name);
                    incorrect = true;
                    break;
                }
                self.items[i].value = Some(arg.clone());
                matched = None;
            } else if let Some(name) = arg.strip_prefix('-') {
                let Some(i) = self.items.iter().position(|it| it.name == name) else {
                    crate::error!("unregistered config item: {arg}");
                    incorrect = true;
                    break;
                };
                if self.items[i].ty == ItemType::Flag {
                    self.items[i].value = Some("true".to_string());
                } else {
                    matched = Some(i);
                }
            } else {
                crate::error!("expecting config name but encountered invalid argument: {arg}");
                incorrect = true;
                break;
            }
        }
        if incorrect {
            for line in self.item_lines() {
                crate::info!("{line}");
            }
            return Err(anyhow!("invalid command line"));
        }
        // ConfigBuilder.isValid: required items given, given paths existing, numbers (given or default) parsing
        let mut ok = true;
        for item in &self.items {
            if item.required && item.value.is_none() {
                crate::error!("missing required config: {}: {}", item.name, item.desc);
                ok = false;
                continue;
            }
            let value = item.value.clone().or_else(|| item.default.clone());
            match item.ty {
                ItemType::Path | ItemType::Paths => {
                    if let Some(v) = &item.value {
                        if !self.valid_path(item, v) {
                            crate::error!("invalid path for config: {} = {}", item.name, v);
                            ok = false;
                        }
                    }
                }
                ItemType::Integer if value.as_deref().and_then(crate::java_text::parse_int).is_none() => {
                    crate::error!("invalid type for config: {}: {} = {}", item.name, item.ty.name(), value.as_deref().unwrap_or("null"));
                    ok = false;
                }
                ItemType::Decimal if value.as_deref().and_then(crate::java_text::parse_double).is_none() => {
                    crate::error!("invalid type for config: {}: {} = {}", item.name, item.ty.name(), value.as_deref().unwrap_or("null"));
                    ok = false;
                }
                _ => {}
            }
        }
        if !ok {
            return Err(anyhow!("invalid command line"));
        }
        // ConfigUtils.setLogLevel
        let debug = self.is_registered("log_debug") && self.has_flag("log_debug");
        let level = if self.is_registered("log_level") { self.value("log_level").unwrap_or_default().to_ascii_uppercase() } else { String::new() };
        if debug || level == "DEBUG" || level == "TRACE" || level == "ALL" {
            crate::log::DEBUG.store(true, std::sync::atomic::Ordering::Relaxed);
        }
        crate::info!("{} version {}", self.app_name, self.version);
        Ok(Parsed::Run)
    }
}

/// ConfigUtils.convertWildcardSamplePath(path, sampleId): `samplePath.replaceAll("\\*", sampleId)`, the
/// sample id read as a Matcher replacement string ('\' escapes, '$0' the match, other '$' references
/// rejected), as Java does. (Added for CUPPA from chord-rs's copy.)
pub fn convert_wildcard_sample_path(sample_path: &str, sample_id: &str) -> Result<String> {
    if !sample_path.contains('*') {
        return Ok(sample_path.to_string());
    }
    let replacement = java_replacement(sample_id, "*")?;
    Ok(sample_path.replace('*', &replacement))
}

/// Matcher.appendReplacement's reading of a replacement string for a pattern without groups whose match
/// is `group0`
pub fn java_replacement(replacement: &str, group0: &str) -> Result<String> {
    let chars: Vec<char> = replacement.chars().collect();
    let mut out = String::new();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c == '\\' {
            i += 1;
            if i == chars.len() {
                return Err(anyhow!("java.lang.IllegalArgumentException: character to be escaped is missing"));
            }
            out.push(chars[i]);
            i += 1;
        } else if c == '$' {
            i += 1;
            if i == chars.len() {
                return Err(anyhow!("java.lang.IllegalArgumentException: Illegal group reference: group index is missing"));
            }
            if chars[i] == '{' {
                let name: String = chars[i + 1..].iter().take_while(|c| **c != '}').collect();
                return Err(anyhow!("java.lang.IllegalArgumentException: No group with name {{{name}}}"));
            }
            let Some(d) = chars[i].to_digit(10) else {
                return Err(anyhow!("java.lang.IllegalArgumentException: Illegal group reference"));
            };
            i += 1;
            if d != 0 {
                return Err(anyhow!("java.lang.IndexOutOfBoundsException: No group {d}"));
            }
            out.push_str(group0);
        } else {
            out.push(c);
            i += 1;
        }
    }
    Ok(out)
}

/// FileWriterUtils.checkAddDirSeparator
pub fn check_add_dir_separator(dir: &str) -> String {
    if dir.ends_with('/') { dir.to_string() } else { format!("{dir}/") }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_like_java() {
        let mut cb = ConfigBuilder::new("Test", "1.0");
        cb.add_str("tumor", false, "Tumor ID");
        cb.add_integer("threads", "threads", 1);
        cb.add_decimal("min_purity", "min", 0.08);
        cb.add_flag("no_charts", "no charts");
        let argv: Vec<String> = ["-tumor", "T", "-threads", "-2", "-no_charts"].iter().map(|s| s.to_string()).collect();
        assert!(matches!(cb.parse(&argv), Ok(Parsed::Run)));
        assert_eq!(cb.value("tumor").as_deref(), Some("T"));
        assert_eq!(cb.integer("threads"), -2);
        assert_eq!(cb.decimal("min_purity"), 0.08);
        assert!(cb.has_flag("no_charts"));
        assert_eq!(cb.value("min_purity").as_deref(), Some("0.08"));
        let mut cb2 = ConfigBuilder::new("Test", "1.0");
        cb2.add_str("tumor", false, "Tumor ID");
        assert!(cb2.parse(&["-tumor".to_string(), "-x".to_string()]).is_err());
    }
}
