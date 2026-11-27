//! FileBasedPlotManager: plots copied into `<output_dir>/plot/` and named relative to the output directory.

use anyhow::{bail, Result};

use crate::compat::java;

pub struct PlotManager {
    output_dir: String,
}

impl PlotManager {
    pub fn new(output_dir: &str) -> PlotManager {
        PlotManager { output_dir: output_dir.to_string() }
    }

    /// `outputDir + File.separator + "plot"` (a doubled separator when the output directory ends in one)
    fn plot_dir(&self) -> String {
        format!("{}/plot", self.output_dir)
    }

    /// createPlotDirectory
    pub fn create_plot_directory(&self) -> Result<()> {
        let dir = self.plot_dir();
        let p = std::path::Path::new(&dir);
        if p.exists() {
            if !p.is_dir() {
                bail!("java.lang.IllegalStateException: Plot directory of [{dir}] is not a directory. Please check configured plot directory.");
            }
        } else {
            std::fs::create_dir_all(p).map_err(|e| anyhow::anyhow!("java.io.IOException: Unable to create plot directory: {dir} ({e})"))?;
            hmf_common::debug!("Created plot directory '{dir}'");
        }
        Ok(())
    }

    /// processPlotFile: None for no plot or a missing file (with a warning); the plot is copied unless the
    /// target exists, and named relative to the output directory: "/plot/<file>" in Java (OR2), "plot/<file>"
    pub fn process_plot_file(&self, source: Option<&str>) -> Result<Option<String>> {
        let Some(source) = source else { return Ok(None) };
        if !std::path::Path::new(source).exists() {
            hmf_common::warn!("Missing source plot path: {source}");
            return Ok(None);
        }
        let name = &source[source.rfind('/').map(|i| i + 1).unwrap_or(0)..];
        let target = format!("{}/{name}", self.plot_dir());
        if !std::path::Path::new(&target).exists() {
            std::fs::copy(source, &target).map_err(|e| anyhow::anyhow!("copying {source} to {target}: {e}"))?;
        }
        Ok(Some(self.relative_path(&target)?))
    }

    /// relativePath(target, outputDir)
    fn relative_path(&self, target: &str) -> Result<String> {
        let root = &self.output_dir;
        let remove = if root.ends_with('/') { root.clone() } else { format!("{root}/") };
        let Some(i) = target.find(&remove) else {
            bail!("java.lang.IllegalStateException: Cannot make relative path of '{target}' based on '{root}'");
        };
        let rel = &target[i + remove.len()..];
        Ok(if java("OR2") { rel.to_string() } else { rel.trim_start_matches('/').to_string() })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn relative_paths() {
        // FileBasedPlotManagerTest
        let m = PlotManager::new("/out/");
        crate::compat::with_java(true, || {
            assert_eq!(m.relative_path("/out//plot/a.png").unwrap(), "/plot/a.png");
            assert_eq!(m.relative_path("/out/plot/a.png").unwrap(), "plot/a.png");
            assert!(m.relative_path("/elsewhere/a.png").is_err());
        });
        crate::compat::with_java(false, || {
            assert_eq!(m.relative_path("/out//plot/a.png").unwrap(), "plot/a.png");
        });
    }
}
