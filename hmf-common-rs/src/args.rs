//! Command line arguments in the style of the Java ConfigBuilder: `-name value` and `-flag`.

use anyhow::{bail, Result};
use std::collections::HashMap;

pub struct Args {
    values: HashMap<String, String>,
    flags: Vec<String>,
}

impl Args {
    /// `flags` lists the options that take no value; every other option takes one.
    pub fn parse(args: &[String], flags: &[&str], items: &[&str]) -> Result<Args> {
        let mut values = HashMap::new();
        let mut found = Vec::new();
        let mut i = 0;
        while i < args.len() {
            let a = &args[i];
            let Some(name) = a.strip_prefix('-') else { bail!("unexpected argument: {a}") };
            if flags.contains(&name) {
                found.push(name.to_string());
                i += 1;
            } else if items.contains(&name) {
                if i + 1 >= args.len() {
                    bail!("missing value for -{name}");
                }
                values.insert(name.to_string(), args[i + 1].clone());
                i += 2;
            } else {
                bail!("unknown option: -{name}");
            }
        }
        Ok(Args { values, flags: found })
    }
    pub fn has_flag(&self, name: &str) -> bool {
        self.flags.iter().any(|f| f == name)
    }
    pub fn has_value(&self, name: &str) -> bool {
        self.values.contains_key(name)
    }
    pub fn value(&self, name: &str) -> Option<&str> {
        self.values.get(name).map(|s| s.as_str())
    }
    pub fn value_or<'a>(&'a self, name: &str, default: &'a str) -> &'a str {
        self.value(name).unwrap_or(default)
    }
    pub fn integer(&self, name: &str, default: i32) -> Result<i32> {
        match self.value(name) {
            Some(v) => Ok(v.parse()?),
            None => Ok(default),
        }
    }
    pub fn decimal(&self, name: &str, default: f64) -> Result<f64> {
        match self.value(name) {
            Some(v) => Ok(v.parse()?),
            None => Ok(default),
        }
    }
}
