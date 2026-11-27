//! OrangeConfig: the options of ORANGE 5.1 (OrangeConfig.registerConfig) and how the tool directories are
//! resolved (PathResolver, PipelineToolDirectories, PathUtil).

use anyhow::{anyhow, bail, Result};
use hmf_common::config_builder::{ConfigBuilder, ItemType};

pub const APP_NAME: &str = "Orange";
pub const VERSION: &str = "5.1";

/// (display_sample_id and add_disclaimer are used only by the PDF, which the port does not write)
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct OrangeConfig {
    /// ExperimentType name: WHOLE_GENOME or TARGETED
    pub run_type: String,
    pub tumor_id: String,
    pub reference_id: Option<String>,
    pub rna_sample_id: Option<String>,
    pub display_sample_id: String,
    pub panel_name: Option<String>,
    pub ref_genome_v37: bool,
    /// SequencingType name
    pub sequencing_type: String,
    /// in the order Java's HashSet gives them
    pub primary_tumor_doids: Vec<String>,
    pub primary_tumor_location: Option<String>,
    /// yyyy-MM-dd
    pub sampling_date: String,
    pub output_dir: String,
    pub output_id: Option<String>,
    pub doid_json_file: Option<String>,
    pub pipeline_version_file: Option<String>,
    pub purple_data_directory: String,
    pub purple_rna_append_vcf: Option<String>,
    pub purple_plot_directory: String,
    pub qsee_directory: Option<String>,
    pub linx_somatic_data_directory: String,
    pub linx_germline_data_directory: Option<String>,
    pub linx_plot_directory: Option<String>,
    pub sage_plot_directory: Option<String>,
    pub lilac_dir: Option<String>,
    pub chord_dir: Option<String>,
    pub cuppa_dir: Option<String>,
    pub peach_dir: Option<String>,
    pub sigs_dir: Option<String>,
    pub virus_dir: Option<String>,
    pub isofox_dir: Option<String>,
    pub add_disclaimer: bool,
    /// the port's -no_pdf (the PDF is never written; accepted for clarity)
    pub no_pdf: bool,
}

impl OrangeConfig {
    pub fn has_reference(&self) -> bool {
        self.reference_id.is_some()
    }
    pub fn has_rna(&self) -> bool {
        self.rna_sample_id.is_some() || self.isofox_dir.is_some()
    }
    pub fn is_targeted(&self) -> bool {
        self.run_type == "TARGETED"
    }
}

const TUMOR_DESC: &str = "Tumor sample ID";
const REFERENCE_DESC: &str = "Reference sample ID";

fn tool_dir_desc(tool: &str) -> String {
    format!("Path to {tool} pipeline files")
}

/// OrangeConfig.registerConfig (and the port's options)
pub fn register(cb: &mut ConfigBuilder) {
    cb.add_str("experiment_type", true, "The type of the experiment, one of WGS or PANEL");
    cb.add_str("tumor", true, TUMOR_DESC);
    cb.add_str("reference", false, REFERENCE_DESC);
    cb.add_str("primary_tumor_location", false, "Primary tumor location displayed in the report");
    cb.add_str("primary_tumor_doids", false, "A semicolon-separated list of DOIDs representing the primary tumor of patient");
    cb.add_str("sampling_date", false, "Optional, if provided represents the sampling date in YYMMDD format");
    cb.add_str_default("ref_genome_version", false, "Ref genome version, 37 or 38", "V37");
    cb.add_str_default("sequencing_type", false, "Sequencing types: ILLUMINA, ULTIMA, SBX", "ILLUMINA");
    cb.add_str("output_dir", false, "Output directory");
    cb.add_str("output_id", false, "Output file suffix");
    cb.add_path("doid_json", false, "Path to JSON file containing the full DOID tree");
    cb.add_path("pipeline_version_file", false, "Path towards the pipeline version file.");
    cb.add_path("pipeline_sample_root_dir", false, "Pipeline sample root directory");
    cb.add_path("sample_data_dir", false, "Sample data directory for all inputs");
    for (name, tool) in [
        ("linx_germline_dir", "linx_germline"),
        ("virus_dir", "virus"),
        ("chord_dir", "chord"),
        ("cuppa_dir", "cuppa"),
        ("peach_dir", "peach"),
        ("sigs_dir", "sigs"),
        ("sage_dir", "sage"),
        ("sage_germline_dir", "sage_germline"),
        ("purple_dir", "purple"),
        ("purple_plot_dir", "purple_plot"),
        ("linx_dir", "linx"),
        ("linx_plot_dir", "linx_plot"),
        ("sage_plot_dir", "sage_plot"),
        ("lilac_dir", "lilac"),
        ("qsee_dir", "qsee"),
    ] {
        cb.add_path(name, false, &tool_dir_desc(tool));
    }
    cb.add_str(
        "pipeline_format",
        false,
        "Assumed directory structure for tool directories. Possible values: OA_V3_1, OA_V3_0, OA_V2_3, OA_V2_2, OA_V2_0, PIP5_V6_0, DB_V6_0. Default: OA_V3_1",
    );
    cb.add_path("pipeline_format_file", false, "File describing expected tool directory structure.");
    cb.add_str("display_sample", false, "Sample ID for display instead of Tumor ID for demo purposes");
    cb.add_str("panel_name", false, "Panel name if known");
    cb.add_str("rna_sample", false, "RNA sample ID");
    cb.add_path("isofox_dir", false, &tool_dir_desc("isofox"));
    cb.add_path("purple_rna_append_vcf", false, "Purple somatic VCF appended with RNA");
    cb.add_flag("add_disclaimer", "Prints a disclaimer on each page");
    cb.add_logging_options();
    cb.add(ItemType::Flag, "no_pdf", false, "Port option: do not write the PDF (the port never writes it)", None);
}

/// FileWriterUtils.checkAddDirSeparator (null and "" unchanged)
pub fn check_add_dir_separator(dir: Option<&str>) -> Option<String> {
    dir.map(|d| if d.is_empty() || d.ends_with('/') { d.to_string() } else { format!("{d}/") })
}

fn exists(path: &str) -> bool {
    std::path::Path::new(path).exists()
}

/// The tool directories of a pipeline layout (PipelineToolDirectories), those ORANGE uses
#[derive(Clone, Debug, Default)]
struct ToolDirectories {
    chord: String,
    cuppa: String,
    isofox: String,
    lilac: String,
    linx_germline: String,
    linx_somatic: String,
    peach: String,
    purple: String,
    sigs: String,
    virus_interpreter: String,
    qsee: String,
}

fn layout(format: &str) -> Result<ToolDirectories> {
    let d = |chord: &str, linx_germline: &str, linx_somatic: &str, virus: &str, qsee: &str| ToolDirectories {
        chord: chord.into(),
        cuppa: "cuppa".into(),
        isofox: "isofox".into(),
        lilac: "lilac".into(),
        linx_germline: linx_germline.into(),
        linx_somatic: linx_somatic.into(),
        peach: "peach".into(),
        purple: "purple".into(),
        sigs: "sigs".into(),
        virus_interpreter: virus.into(),
        qsee: qsee.into(),
    };
    Ok(match format {
        "OA_V3_1" | "OA_V3_0" | "OA_V2_3" | "OA_V2_2" | "DB_V6_0" => {
            d("chord", "linx/germline_annotations", "linx/somatic_annotations", "virusinterpreter", "qsee")
        }
        "OA_V2_0" => d("chord", "linx/germline_annotations", "linx/somatic_annotations", "virusinterpreter", ""),
        "PIP5_V6_0" => d("chord", "linx_germline", "linx", "virusintrprtr", "qsee"),
        other => bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.pipeline.PipelineOutputStructure.{other}"),
    })
}

/// PipelineToolDirectoriesFile.read: `name<TAB>dir` lines, "" for an absent name
fn layout_from_file(path: &str) -> Result<ToolDirectories> {
    let lines = hmf_common::java_text::read_lines(path)
        .map_err(|e| anyhow!("java.lang.IllegalArgumentException: Could not load tool subdirectory config file: {path} ({e})"))?;
    let get = |field: &str| -> Result<String> {
        for line in &lines {
            let v: Vec<&str> = line.split('\t').collect();
            if v[0] == field {
                return v.get(1).map(|s| s.to_string()).ok_or_else(|| anyhow!("java.lang.ArrayIndexOutOfBoundsException: {path}"));
            }
        }
        Ok(String::new())
    };
    Ok(ToolDirectories {
        chord: get("chordDir")?,
        cuppa: get("cuppaDir")?,
        isofox: get("isofoxDir")?,
        lilac: get("lilacDir")?,
        linx_germline: get("linxGermlineDir")?,
        linx_somatic: get("linxSomaticDir")?,
        peach: get("peachDir")?,
        purple: get("purpleDir")?,
        sigs: get("sigsDir")?,
        virus_interpreter: get("virusInterpreterDir")?,
        qsee: get("qSeeDir")?,
    })
}

/// ConfigUtils.convertWildcardSamplePath ('$' the normal, when there is one, then '*' the tumor)
fn convert_wildcards(path: &str, tumor: &str, normal: Option<&str>) -> String {
    match normal {
        Some(n) => path.replace('$', n).replace('*', tumor),
        None => path.replace('*', tumor),
    }
}

struct PathResolver<'a> {
    cb: &'a ConfigBuilder,
    root: Option<String>,
    sample_data_dir: Option<String>,
}

impl PathResolver<'_> {
    fn optional_tool_dir(&self, key: &str, default_dir: &str) -> Option<String> {
        if self.cb.has_value(key) {
            return check_add_dir_separator(self.cb.value(key).as_deref());
        }
        if let Some(root) = &self.root {
            return check_add_dir_separator(Some(&format!("{root}{default_dir}")));
        }
        self.sample_data_dir.clone()
    }

    fn mandatory_tool_dir(&self, key: &str, default_dir: &str) -> Result<String> {
        let dir = self.optional_tool_dir(key, default_dir).ok_or_else(|| {
            anyhow!("java.lang.IllegalArgumentException: Failed to determine tool directory for configuration [{key}/{default_dir}].")
        })?;
        mandatory_path(&dir)?;
        Ok(check_add_dir_separator(Some(&dir)).unwrap())
    }

    fn optional_tool_plots_dir(&self, key: &str, default_dir: &str) -> Option<String> {
        if self.cb.has_value(key) {
            return check_add_dir_separator(self.cb.value(key).as_deref());
        }
        self.root.as_ref().map(|root| format!("{root}{default_dir}/plot/"))
    }

    fn mandatory_tool_plots_dir(&self, key: &str, default_dir: &str) -> Result<String> {
        let dir = self.optional_tool_plots_dir(key, default_dir).ok_or_else(|| {
            anyhow!("java.lang.IllegalArgumentException: Failed to determine plot directory for configuration [{key}/{default_dir}].")
        })?;
        mandatory_path(&dir)?;
        Ok(check_add_dir_separator(Some(&dir)).unwrap())
    }
}

/// PathUtil.mandatoryPath
fn mandatory_path(path: &str) -> Result<()> {
    if !exists(path) {
        bail!("java.lang.IllegalArgumentException: Path [{path}] did not exist and is mandatory.");
    }
    Ok(())
}

/// PathUtil.optionalPath
fn optional_path(path: Option<String>) -> Option<String> {
    path.filter(|p| exists(p))
}

/// The local date as LocalDate.now() gives it (the system time zone)
pub fn local_date_today() -> (i32, u32, u32) {
    unsafe extern "C" {
        fn localtime_r(t: *const i64, tm: *mut [i32; 16]) -> *mut [i32; 16];
    }
    let secs = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs() as i64;
    let mut raw = [0i32; 16];
    unsafe { localtime_r(&secs, &mut raw) };
    (raw[5] + 1900, (raw[4] + 1) as u32, raw[3] as u32)
}

fn is_leap(y: i32) -> bool {
    (y % 4 == 0 && y % 100 != 0) || y % 400 == 0
}

fn month_length(y: i32, m: u32) -> u32 {
    match m {
        2 => {
            if is_leap(y) {
                29
            } else {
                28
            }
        }
        4 | 6 | 9 | 11 => 30,
        _ => 31,
    }
}

/// LocalDate.parse(s, ofPattern("yyMMdd", ENGLISH)): two-digit years 2000-2099, ResolverStyle.SMART (a
/// day of 29 to 31 past the month's end becomes its last day); None where Java throws DateTimeParseException
pub fn parse_sampling_date(s: &str) -> Option<(i32, u32, u32)> {
    if s.len() != 6 || !s.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let y = 2000 + s[0..2].parse::<i32>().ok()?;
    let m = s[2..4].parse::<u32>().ok()?;
    let d = s[4..6].parse::<u32>().ok()?;
    if !(1..=12).contains(&m) || !(1..=31).contains(&d) {
        return None;
    }
    Some((y, m, d.min(month_length(y, m))))
}

pub fn format_date((y, m, d): (i32, u32, u32)) -> String {
    format!("{y:04}-{m:02}-{d:02}")
}

/// Sets.newHashSet(...) of strings, in Java's iteration order
fn java_hash_set_order(values: Vec<String>) -> Vec<String> {
    let mut set = hmf_common::java_hash::JavaHashSet::new();
    for v in values {
        set.insert(v);
    }
    set.to_vec()
}

impl OrangeConfig {
    /// new OrangeConfig(configBuilder)
    pub fn from_builder(cb: &ConfigBuilder) -> Result<OrangeConfig> {
        let tumor_id = cb.value("tumor").unwrap();
        let reference_id = cb.value("reference");

        let (primary_tumor_doids, doid_json_file) = if cb.has_value("primary_tumor_doids") && cb.has_value("doid_json") {
            let values: Vec<String> = cb.value("primary_tumor_doids").unwrap().split(';').map(String::from).collect();
            (java_hash_set_order(values), cb.value("doid_json"))
        } else {
            (Vec::new(), None)
        };
        let primary_tumor_location = cb.value("primary_tumor_location");

        let ref_genome_v37 = match cb.value("ref_genome_version").as_deref() {
            None | Some("V37") | Some("37") | Some("HG37") => true,
            Some("V38") | Some("38") | Some("HG38") => false,
            Some(v) => bail!("java.lang.IllegalArgumentException: Cannot resolve ref genome version: {v}"),
        };
        let sequencing_type = match cb.value("sequencing_type") {
            None => "ILLUMINA".to_string(),
            Some(v) if ["ILLUMINA", "ULTIMA", "SBX"].contains(&v.as_str()) => v,
            Some(v) => bail!("java.lang.IllegalArgumentException: No enum constant com.hartwig.hmftools.common.sequencing.SequencingType.{v}"),
        };

        let root = check_add_dir_separator(cb.value("pipeline_sample_root_dir").as_deref());
        let mut pipeline_version_file = cb.value("pipeline_version_file");
        if pipeline_version_file.is_none() {
            if let Some(r) = &root {
                let test = format!("{r}/orange_pipeline.version.txt");
                if exists(&test) {
                    pipeline_version_file = Some(test);
                }
            }
        }

        let output_dir = if root.is_some() && !cb.has_value("output_dir") {
            format!("{}orange/", root.as_ref().unwrap())
        } else {
            check_add_dir_separator(cb.value("output_dir").as_deref())
                .ok_or_else(|| anyhow!("java.lang.NullPointerException: no output directory (-output_dir)"))?
        };
        if !exists(&output_dir) {
            let _ = std::fs::create_dir_all(&output_dir);
        }
        let output_id = cb.value("output_id");

        let resolver = PathResolver {
            cb,
            root: root.clone(),
            sample_data_dir: check_add_dir_separator(cb.value("sample_data_dir").as_deref()),
        };

        // resolveDefaultPipelineDirectories: always resolved (an invalid -pipeline_format fails even when
        // every directory is given)
        let mut dirs = if cb.has_value("pipeline_format_file") {
            layout_from_file(&cb.value("pipeline_format_file").unwrap())?
        } else {
            layout(cb.value("pipeline_format").as_deref().unwrap_or("OA_V3_1"))?
        };
        {
            let n = reference_id.as_deref();
            for d in [
                &mut dirs.chord,
                &mut dirs.cuppa,
                &mut dirs.isofox,
                &mut dirs.lilac,
                &mut dirs.linx_germline,
                &mut dirs.linx_somatic,
                &mut dirs.peach,
                &mut dirs.purple,
                &mut dirs.sigs,
                &mut dirs.virus_interpreter,
                &mut dirs.qsee,
            ] {
                *d = convert_wildcards(d, &tumor_id, n);
            }
        }

        let purple_data_directory = resolver.mandatory_tool_dir("purple_dir", &dirs.purple)?;
        let purple_plot_directory = resolver.mandatory_tool_plots_dir("purple_plot_dir", &dirs.purple)?;
        let linx_somatic_data_directory = resolver.mandatory_tool_dir("linx_dir", &dirs.linx_somatic)?;
        let linx_plot_directory = optional_path(resolver.optional_tool_plots_dir("linx_plot_dir", &dirs.linx_somatic));
        // the SAGE plot directory's default is under LINX's somatic directory, as in Java
        let sage_plot_directory = optional_path(resolver.optional_tool_plots_dir("sage_plot_dir", &dirs.linx_somatic));
        let qsee_directory = resolver.optional_tool_dir("qsee_dir", &dirs.qsee);
        let linx_germline_data_directory = if reference_id.is_some() {
            Some(resolver.mandatory_tool_dir("linx_germline_dir", &dirs.linx_germline)?)
        } else {
            None
        };
        let lilac_dir = resolver.optional_tool_dir("lilac_dir", &dirs.lilac);
        let chord_dir = resolver.optional_tool_dir("chord_dir", &dirs.chord);
        let cuppa_dir = resolver.optional_tool_dir("cuppa_dir", &dirs.cuppa);
        let peach_dir = resolver.optional_tool_dir("peach_dir", &dirs.peach);
        let sigs_dir = resolver.optional_tool_dir("sigs_dir", &dirs.sigs);
        let virus_dir = resolver.optional_tool_dir("virus_dir", &dirs.virus_interpreter);

        let rna_sample_id = cb.value("rna_sample");
        let purple_rna_append_vcf = cb.value("purple_rna_append_vcf");
        let isofox_dir = if cb.has_value("isofox_dir") {
            cb.value("isofox_dir")
        } else {
            resolver.optional_tool_dir("isofox_dir", &dirs.isofox).filter(|d| exists(d))
        };
        if isofox_dir.is_some() || rna_sample_id.is_some() {
            hmf_common::debug!("RNA sample({:?}) Isofox results({:?})", rna_sample_id, isofox_dir);
        }

        let display_sample_id = cb.value("display_sample").unwrap_or_else(|| tumor_id.clone());
        let panel_name = cb.value("panel_name");
        let add_disclaimer = cb.has_flag("add_disclaimer");

        let sampling_date = match cb.value("sampling_date") {
            Some(s) => match parse_sampling_date(&s) {
                Some(d) => format_date(d),
                None => {
                    hmf_common::warn!("Could not parse configured sampling date '{s}'. Expected format is 'yyMMdd'");
                    format_date(local_date_today())
                }
            },
            None => format_date(local_date_today()),
        };

        let run_type = match cb.value("experiment_type").as_deref() {
            Some("WGS") => "WHOLE_GENOME",
            Some("PANEL") => "TARGETED",
            _ => bail!("java.lang.IllegalArgumentException: Invalid experiment type, must be one of WGS, PANEL"),
        }
        .to_string();
        hmf_common::info!("experiment type has been resolved to '{run_type}'");

        Ok(OrangeConfig {
            run_type,
            tumor_id,
            reference_id,
            rna_sample_id,
            display_sample_id,
            panel_name,
            ref_genome_v37,
            sequencing_type,
            primary_tumor_doids,
            primary_tumor_location,
            sampling_date,
            output_dir,
            output_id,
            doid_json_file,
            pipeline_version_file,
            purple_data_directory,
            purple_rna_append_vcf,
            purple_plot_directory,
            qsee_directory,
            linx_somatic_data_directory,
            linx_germline_data_directory,
            linx_plot_directory,
            sage_plot_directory,
            lilac_dir,
            chord_dir,
            cuppa_dir,
            peach_dir,
            sigs_dir,
            virus_dir,
            isofox_dir,
            add_disclaimer,
            no_pdf: cb.has_flag("no_pdf"),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sampling_dates() {
        assert_eq!(parse_sampling_date("230131").map(format_date).as_deref(), Some("2023-01-31"));
        // SMART resolution: the 30th of February is its last day
        assert_eq!(parse_sampling_date("240230").map(format_date).as_deref(), Some("2024-02-29"));
        assert_eq!(parse_sampling_date("230230").map(format_date).as_deref(), Some("2023-02-28"));
        assert!(parse_sampling_date("231301").is_none());
        assert!(parse_sampling_date("230132").is_none());
        assert!(parse_sampling_date("2023-01-31").is_none());
        assert!(parse_sampling_date("23013").is_none());
    }

    #[test]
    fn wildcards() {
        assert_eq!(convert_wildcards("bamtools/$_x/*", "T", Some("N")), "bamtools/N_x/T");
        assert_eq!(convert_wildcards("bamtools/$_x/*", "T", None), "bamtools/$_x/T");
    }
}
