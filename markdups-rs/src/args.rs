//! MarkDuplicates' arguments as GATK 4.6.1.0 (Barclay) parses them, and the command line it records in
//! the @PG line and the metrics file: the arguments given, in the order the fields are declared in
//! Picard's classes, then the defaults of the others in the same order (null values and empty lists
//! left out).

use anyhow::{Result, bail};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind {
    Bool,
    Int,
    Long,
    Double,
    Str,
    /// java.io.File: printed in its normalised form
    File,
    Enum(&'static [&'static str]),
    StrList,
    FileList,
}

/// What the port does with an argument.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Support {
    /// implemented
    Yes,
    /// accepted; it only tunes Picard's memory use or logging and changes no output
    NoEffect,
    /// accepted only at its default value
    DefaultOnly,
}

pub struct Def {
    pub name: &'static str,
    pub short: Option<&'static str>,
    pub kind: Kind,
    /// the default as Barclay prints it; None for null (or an empty list)
    pub default: Option<&'static str>,
    pub support: Support,
}

const SCORING: &[&str] = &["SUM_OF_BASE_QUALITIES", "TOTAL_MAPPED_REFERENCE_LENGTH", "RANDOM"];
const TAGGING: &[&str] = &["DontTag", "OpticalOnly", "All"];
const SORT_ORDERS: &[&str] = &["unsorted", "queryname", "coordinate", "duplicate", "unknown"];
const FLOW_STRATEGY: &[&str] = &["FLOW_QUALITY_SUM_STRATEGY", "FLOW_END_QUALITY_STRATEGY"];
const LOG_LEVELS: &[&str] = &["ERROR", "WARNING", "INFO", "DEBUG"];
const STRINGENCY: &[&str] = &["STRICT", "LENIENT", "SILENT"];

/// Picard's default read name regex, a marker for the optimised parsing.
pub const DEFAULT_READ_NAME_REGEX: &str = "<optimized capture of last three ':' separated fields as numeric values>";

macro_rules! def {
    ($name:expr, $short:expr, $kind:expr, $default:expr, $support:expr) => {
        Def { name: $name, short: $short, kind: $kind, default: $default, support: $support }
    };
}

/// MarkDuplicates' arguments in declaration order: MarkDuplicates, its flow argument collection,
/// AbstractMarkDuplicatesCommandLineProgram (the PG tag collection first),
/// AbstractOpticalDuplicateFinderCommandLineProgram, CommandLineProgram.
pub static DEFS: &[Def] = &[
    def!("MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP", Some("MAX_SEQS"), Kind::Int, Some("50000"), Support::NoEffect),
    def!("MAX_FILE_HANDLES_FOR_READ_ENDS_MAP", Some("MAX_FILE_HANDLES"), Kind::Int, Some("8000"), Support::NoEffect),
    def!("SORTING_COLLECTION_SIZE_RATIO", None, Kind::Double, Some("0.25"), Support::NoEffect),
    def!("BARCODE_TAG", None, Kind::Str, None, Support::DefaultOnly),
    def!("READ_ONE_BARCODE_TAG", None, Kind::Str, None, Support::DefaultOnly),
    def!("READ_TWO_BARCODE_TAG", None, Kind::Str, None, Support::DefaultOnly),
    def!("TAG_DUPLICATE_SET_MEMBERS", None, Kind::Bool, Some("false"), Support::DefaultOnly),
    def!("REMOVE_SEQUENCING_DUPLICATES", None, Kind::Bool, Some("false"), Support::Yes),
    def!("TAGGING_POLICY", None, Kind::Enum(TAGGING), Some("DontTag"), Support::Yes),
    def!("CLEAR_DT", None, Kind::Bool, Some("true"), Support::Yes),
    def!("DUPLEX_UMI", None, Kind::Bool, Some("false"), Support::DefaultOnly),
    def!("MOLECULAR_IDENTIFIER_TAG", None, Kind::Str, None, Support::DefaultOnly),
    def!("FLOW_MODE", None, Kind::Bool, Some("false"), Support::DefaultOnly),
    def!("FLOW_DUP_STRATEGY", None, Kind::Enum(FLOW_STRATEGY), Some("FLOW_QUALITY_SUM_STRATEGY"), Support::NoEffect),
    def!("FLOW_USE_END_IN_UNPAIRED_READS", None, Kind::Bool, Some("false"), Support::NoEffect),
    def!("FLOW_USE_UNPAIRED_CLIPPED_END", None, Kind::Bool, Some("false"), Support::NoEffect),
    def!("FLOW_UNPAIRED_END_UNCERTAINTY", None, Kind::Int, Some("0"), Support::NoEffect),
    def!("FLOW_UNPAIRED_START_UNCERTAINTY", None, Kind::Int, Some("0"), Support::NoEffect),
    def!("FLOW_SKIP_FIRST_N_FLOWS", None, Kind::Int, Some("0"), Support::NoEffect),
    def!("FLOW_Q_IS_KNOWN_END", None, Kind::Bool, Some("false"), Support::NoEffect),
    def!("FLOW_EFFECTIVE_QUALITY_THRESHOLD", None, Kind::Int, Some("15"), Support::NoEffect),
    def!("ADD_PG_TAG_TO_READS", None, Kind::Bool, Some("true"), Support::Yes),
    def!("INPUT", Some("I"), Kind::StrList, None, Support::Yes),
    def!("OUTPUT", Some("O"), Kind::File, None, Support::Yes),
    def!("METRICS_FILE", Some("M"), Kind::File, None, Support::Yes),
    def!("REMOVE_DUPLICATES", None, Kind::Bool, Some("false"), Support::Yes),
    def!("ASSUME_SORTED", Some("AS"), Kind::Bool, Some("false"), Support::Yes),
    def!("ASSUME_SORT_ORDER", Some("ASO"), Kind::Enum(SORT_ORDERS), None, Support::Yes),
    def!("DUPLICATE_SCORING_STRATEGY", Some("DS"), Kind::Enum(SCORING), Some("SUM_OF_BASE_QUALITIES"), Support::Yes),
    def!("PROGRAM_RECORD_ID", Some("PG"), Kind::Str, Some("MarkDuplicates"), Support::Yes),
    def!("PROGRAM_GROUP_VERSION", Some("PG_VERSION"), Kind::Str, None, Support::Yes),
    def!("PROGRAM_GROUP_COMMAND_LINE", Some("PG_COMMAND"), Kind::Str, None, Support::Yes),
    def!("PROGRAM_GROUP_NAME", Some("PG_NAME"), Kind::Str, Some("MarkDuplicates"), Support::Yes),
    def!("COMMENT", Some("CO"), Kind::StrList, None, Support::Yes),
    def!("READ_NAME_REGEX", None, Kind::Str, Some(DEFAULT_READ_NAME_REGEX), Support::Yes),
    def!("OPTICAL_DUPLICATE_PIXEL_DISTANCE", None, Kind::Int, Some("100"), Support::Yes),
    def!("MAX_OPTICAL_DUPLICATE_SET_SIZE", None, Kind::Long, Some("300000"), Support::Yes),
    def!("TMP_DIR", None, Kind::FileList, None, Support::Yes),
    def!("VERBOSITY", None, Kind::Enum(LOG_LEVELS), Some("INFO"), Support::NoEffect),
    def!("QUIET", None, Kind::Bool, Some("false"), Support::NoEffect),
    def!("VALIDATION_STRINGENCY", None, Kind::Enum(STRINGENCY), Some("STRICT"), Support::Yes),
    def!("COMPRESSION_LEVEL", None, Kind::Int, Some("2"), Support::Yes),
    def!("MAX_RECORDS_IN_RAM", None, Kind::Int, Some("500000"), Support::NoEffect),
    def!("CREATE_INDEX", None, Kind::Bool, Some("false"), Support::Yes),
    def!("CREATE_MD5_FILE", None, Kind::Bool, Some("false"), Support::DefaultOnly),
    def!("REFERENCE_SEQUENCE", Some("R"), Kind::File, None, Support::NoEffect),
    def!("help", Some("h"), Kind::Bool, Some("false"), Support::DefaultOnly),
    def!("version", None, Kind::Bool, Some("false"), Support::DefaultOnly),
    def!("showHidden", None, Kind::Bool, Some("false"), Support::DefaultOnly),
    def!("USE_JDK_DEFLATER", Some("use_jdk_deflater"), Kind::Bool, Some("false"), Support::NoEffect),
    def!("USE_JDK_INFLATER", Some("use_jdk_inflater"), Kind::Bool, Some("false"), Support::NoEffect),
];

/// The value of one argument as parsed: None for null.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    Null,
    One(String),
    List(Vec<String>),
}

pub struct Parsed {
    /// per definition: Some(value) if given on the command line
    pub given: Vec<Option<Value>>,
    /// the default COMPRESSION_LEVEL (samjdk.compression_level)
    pub compression_default: i32,
}

pub fn index_of(name: &str) -> usize {
    DEFS.iter().position(|d| d.name == name).unwrap_or_else(|| panic!("no argument {name}"))
}

/// java.io.File's normalised path: repeated separators collapsed, a trailing one dropped.
pub fn normalize_file(p: &str) -> String {
    let mut out = String::with_capacity(p.len());
    let mut last_slash = false;
    for c in p.chars() {
        if c == '/' {
            if last_slash {
                continue;
            }
            last_slash = true;
        } else {
            last_slash = false;
        }
        out.push(c);
    }
    if out.len() > 1 && out.ends_with('/') {
        out.pop();
    }
    out
}

fn canonical(kind: Kind, name: &str, v: &str) -> Result<String> {
    Ok(match kind {
        Kind::Bool => match v.to_ascii_lowercase().as_str() {
            "true" => "true".into(),
            "false" => "false".into(),
            _ => bail!("Argument '{name}' has a bad value: {v}. It should be true or false"),
        },
        Kind::Int => match v.trim().parse::<i32>() {
            Ok(i) => i.to_string(),
            Err(_) => bail!("Argument '{name}' has a bad value: {v}. It should be an integer"),
        },
        Kind::Long => match v.trim().parse::<i64>() {
            Ok(i) => i.to_string(),
            Err(_) => bail!("Argument '{name}' has a bad value: {v}. It should be an integer"),
        },
        Kind::Double => match v.trim().parse::<f64>() {
            Ok(d) => hmf_common::fmt::java_double_to_string(d),
            Err(_) => bail!("Argument '{name}' has a bad value: {v}. It should be a number"),
        },
        Kind::Enum(values) => {
            if !values.contains(&v) {
                bail!("Argument '{name}' has a bad value: {v}. Possible values: {{{}}}", values.join(", "));
            }
            v.to_string()
        }
        Kind::File | Kind::FileList => normalize_file(v),
        Kind::Str | Kind::StrList => v.to_string(),
    })
}

/// Parses MarkDuplicates' arguments (after the tool name). `--NAME value` or `-SHORT value`; a
/// boolean without a value is true; "null" sets an argument to null (and empties a list).
pub fn parse(args: &[String], compression_default: i32) -> Result<Parsed> {
    let mut given: Vec<Option<Value>> = vec![None; DEFS.len()];
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        let name = if let Some(n) = a.strip_prefix("--") {
            n
        } else if let Some(n) = a.strip_prefix('-') {
            n
        } else {
            bail!("Unexpected positional argument: {a}");
        };
        if name.contains('=') {
            bail!("Can't parse option name containing an embedded '=' ({name})");
        }
        // Barclay takes either name after either prefix
        let Some(di) = DEFS.iter().position(|d| d.name == name).or_else(|| DEFS.iter().position(|d| d.short == Some(name))) else {
            bail!("{a} is not a recognized option");
        };
        let def = &DEFS[di];
        i += 1;
        let value: Option<String> = if def.kind == Kind::Bool {
            // a boolean takes the next word as its value unless it is an option
            match args.get(i) {
                Some(v) if !v.starts_with('-') => {
                    i += 1;
                    Some(v.clone())
                }
                _ => Some("true".into()),
            }
        } else {
            let Some(v) = args.get(i) else { bail!("Argument '{}' requires a value", def.name) };
            i += 1;
            Some(v.clone())
        };
        let v = value.unwrap();
        let is_list = matches!(def.kind, Kind::StrList | Kind::FileList);
        if v == "null" {
            given[di] = Some(if is_list { Value::List(Vec::new()) } else { Value::Null });
            continue;
        }
        let c = canonical(def.kind, def.name, &v)?;
        if is_list {
            match &mut given[di] {
                Some(Value::List(l)) => l.push(c),
                g => *g = Some(Value::List(vec![c])),
            }
        } else {
            if given[di].is_some() {
                bail!("Argument '{}' cannot be specified more than once", def.name);
            }
            given[di] = Some(Value::One(c));
        }
    }
    for (req, msg) in [("INPUT", "INPUT"), ("OUTPUT", "OUTPUT"), ("METRICS_FILE", "METRICS_FILE")] {
        let g = &given[index_of(req)];
        let missing = match g {
            None | Some(Value::Null) => true,
            Some(Value::List(l)) => l.is_empty(),
            _ => false,
        };
        if missing {
            bail!("Argument '{msg}' is required");
        }
    }
    if given[index_of("ASSUME_SORTED")].is_some() && given[index_of("ASSUME_SORT_ORDER")].is_some() {
        bail!("Argument 'ASSUME_SORTED' cannot be used in conjunction with argument(s) ASSUME_SORT_ORDER");
    }
    Ok(Parsed { given, compression_default })
}

impl Parsed {
    fn default_of(&self, i: usize) -> Option<String> {
        if DEFS[i].name == "COMPRESSION_LEVEL" {
            return Some(self.compression_default.to_string());
        }
        DEFS[i].default.map(str::to_string)
    }
    /// The value of a single-valued argument (given or default); None for null.
    pub fn get(&self, name: &str) -> Option<String> {
        let i = index_of(name);
        match &self.given[i] {
            Some(Value::One(v)) => Some(v.clone()),
            Some(Value::Null) => None,
            Some(Value::List(_)) => unreachable!(),
            None => self.default_of(i),
        }
    }
    pub fn list(&self, name: &str) -> Vec<String> {
        match &self.given[index_of(name)] {
            Some(Value::List(l)) => l.clone(),
            _ => Vec::new(),
        }
    }
    pub fn flag(&self, name: &str) -> bool {
        self.get(name).as_deref() == Some("true")
    }
    pub fn int(&self, name: &str) -> i64 {
        self.get(name).and_then(|v| v.parse().ok()).unwrap_or(0)
    }
    pub fn is_given(&self, name: &str) -> bool {
        self.given[index_of(name)].is_some()
    }

    /// Barclay's command line: "MarkDuplicates", the given arguments, then the defaults. A given null
    /// prints as an empty word.
    pub fn command_line(&self) -> String {
        let mut words: Vec<String> = vec!["MarkDuplicates".into()];
        for (i, d) in DEFS.iter().enumerate() {
            match &self.given[i] {
                Some(Value::One(v)) => words.push(format!("--{} {v}", d.name)),
                Some(Value::List(l)) => {
                    for v in l {
                        words.push(format!("--{} {v}", d.name));
                    }
                }
                Some(Value::Null) => words.push(String::new()),
                None => {}
            }
        }
        for (i, d) in DEFS.iter().enumerate() {
            if self.given[i].is_some() {
                continue;
            }
            if let Some(v) = self.default_of(i) {
                words.push(format!("--{} {v}", d.name));
            }
        }
        words.join(" ")
    }

    /// Fails for arguments the port does not implement, when set to other than their defaults.
    pub fn check_supported(&self) -> Result<()> {
        for (i, d) in DEFS.iter().enumerate() {
            if d.support != Support::DefaultOnly {
                continue;
            }
            let Some(g) = &self.given[i] else { continue };
            let at_default = match g {
                Value::Null => d.default.is_none(),
                Value::One(v) => d.default == Some(v.as_str()),
                Value::List(l) => l.is_empty(),
            };
            if !at_default {
                bail!("--{} is not supported by markdups-rs (only its default value is)", d.name);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(v: &[&str]) -> Vec<String> {
        v.iter().map(|x| x.to_string()).collect()
    }

    #[test]
    fn pipeline_command_line() {
        let p = parse(&s(&["--INPUT", "a.bam", "--OUTPUT", "out/a.md.bam", "--METRICS_FILE", "out/a.md.metrics", "--TMP_DIR", ".", "--CREATE_INDEX"]), 2).unwrap();
        assert_eq!(
            p.command_line(),
            "MarkDuplicates --INPUT a.bam --OUTPUT out/a.md.bam --METRICS_FILE out/a.md.metrics --TMP_DIR . --CREATE_INDEX true --MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP 50000 --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 8000 --SORTING_COLLECTION_SIZE_RATIO 0.25 --TAG_DUPLICATE_SET_MEMBERS false --REMOVE_SEQUENCING_DUPLICATES false --TAGGING_POLICY DontTag --CLEAR_DT true --DUPLEX_UMI false --FLOW_MODE false --FLOW_DUP_STRATEGY FLOW_QUALITY_SUM_STRATEGY --FLOW_USE_END_IN_UNPAIRED_READS false --FLOW_USE_UNPAIRED_CLIPPED_END false --FLOW_UNPAIRED_END_UNCERTAINTY 0 --FLOW_UNPAIRED_START_UNCERTAINTY 0 --FLOW_SKIP_FIRST_N_FLOWS 0 --FLOW_Q_IS_KNOWN_END false --FLOW_EFFECTIVE_QUALITY_THRESHOLD 15 --ADD_PG_TAG_TO_READS true --REMOVE_DUPLICATES false --ASSUME_SORTED false --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES --PROGRAM_RECORD_ID MarkDuplicates --PROGRAM_GROUP_NAME MarkDuplicates --READ_NAME_REGEX <optimized capture of last three ':' separated fields as numeric values> --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 --MAX_OPTICAL_DUPLICATE_SET_SIZE 300000 --VERBOSITY INFO --QUIET false --VALIDATION_STRINGENCY STRICT --COMPRESSION_LEVEL 2 --MAX_RECORDS_IN_RAM 500000 --CREATE_MD5_FILE false --help false --version false --showHidden false --USE_JDK_DEFLATER false --USE_JDK_INFLATER false"
        );
    }

    /// GATK 4.6.1.0 on `-I a.bam --INPUT b.bam -O out3//x.md.bam -M ./out3/x.metrics --READ_NAME_REGEX null
    /// --COMMENT "hello world" --CO two --OPTICAL_DUPLICATE_PIXEL_DISTANCE 0250 --SORTING_COLLECTION_SIZE_RATIO 1e-1
    /// --CREATE_INDEX true --TMP_DIR out3/`
    #[test]
    fn given_order_null_and_normalising() {
        let p = parse(
            &s(&[
                "-I", "a.bam", "--INPUT", "b.bam", "-O", "out3//x.md.bam", "-M", "./out3/x.metrics", "--READ_NAME_REGEX", "null", "--COMMENT", "hello world", "--CO",
                "two", "--OPTICAL_DUPLICATE_PIXEL_DISTANCE", "0250", "--SORTING_COLLECTION_SIZE_RATIO", "1e-1", "--CREATE_INDEX", "true", "--TMP_DIR", "out3/",
            ]),
            2,
        )
        .unwrap();
        let cl = p.command_line();
        assert!(cl.starts_with("MarkDuplicates --SORTING_COLLECTION_SIZE_RATIO 0.1 --INPUT a.bam --INPUT b.bam --OUTPUT out3/x.md.bam --METRICS_FILE ./out3/x.metrics --COMMENT hello world --COMMENT two  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 250 --TMP_DIR out3 --CREATE_INDEX true --MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP 50000"), "{cl}");
        assert!(cl.contains("--PROGRAM_GROUP_NAME MarkDuplicates --MAX_OPTICAL_DUPLICATE_SET_SIZE 300000"));
        assert_eq!(p.get("READ_NAME_REGEX"), None);
    }

    #[test]
    fn errors() {
        assert!(parse(&s(&["--INPUT=a.bam"]), 2).is_err());
        assert!(parse(&s(&["--OUTPUT", "x", "--METRICS_FILE", "m"]), 2).is_err());
        assert!(parse(&s(&["-I", "a", "-O", "x", "-M", "m", "--NOPE", "1"]), 2).is_err());
        let p = parse(&s(&["-I", "a", "-O", "x", "-M", "m", "--BARCODE_TAG", "BC"]), 2).unwrap();
        assert!(p.check_supported().is_err());
    }
}
