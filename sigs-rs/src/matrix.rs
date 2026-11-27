//! hmf-common utils.Matrix and MatrixFile: comma-separated matrices with a header of column names.

use anyhow::{anyhow, bail, Result};

/// A Rows x Cols matrix of doubles (row-major rows)
#[derive(Clone, Debug, PartialEq)]
pub struct Matrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<Vec<f64>>,
}

impl Matrix {
    pub fn new(rows: usize, cols: usize) -> Matrix {
        Matrix { rows, cols, data: vec![vec![0.0; cols]; rows] }
    }

    /// Matrix.get: 0 (and an error logged) for a row or column beyond the matrix
    pub fn get(&self, row: usize, col: usize) -> f64 {
        if row >= self.rows || col >= self.cols {
            hmf_common::error!("incorrect row or column");
            return 0.0;
        }
        self.data[row][col]
    }

    pub fn col(&self, c: usize) -> Vec<f64> {
        self.data.iter().map(|r| r[c]).collect()
    }
}

/// MatrixFile.loadMatrixDataFile(List<String> fileData, columnNames, ignoreFields, transpose).
///
/// `column_names` is the caller's list, to which the header's names are added (and from which the first
/// occurrence of each ignored name is removed), as Java does; the ignored columns are found by their index
/// in that list, so a list that is not empty on entry shifts them (SI4). Rows with another number of
/// fields are left as zeros. Err for "empty matrix data" (Java returns null) and for a value Double.parseDouble
/// rejects (NumberFormatException).
pub fn load_matrix_lines(
    lines: &[String],
    column_names: &mut Vec<String>,
    ignore_fields: Option<&[&str]>,
    transpose: bool,
) -> Result<Matrix> {
    if lines.len() <= 1 {
        hmf_common::error!("empty matrix data");
        bail!("empty matrix data");
    }
    let header = &lines[0];
    column_names.extend(hmf_common::java_text::split_all(header, ',').into_iter().map(String::from));

    let mut ignore_cols: Vec<usize> = Vec::new();
    if let Some(ignore) = ignore_fields {
        for (i, name) in column_names.iter().enumerate() {
            if ignore.contains(&name.as_str()) {
                ignore_cols.push(i);
                if ignore_cols.len() == ignore.len() {
                    break;
                }
            }
        }
        for name in ignore {
            if let Some(p) = column_names.iter().position(|c| c == name) {
                column_names.remove(p);
            }
        }
    }

    let col_count = column_names.len();
    let row_count = lines.len() - 1;
    let mut invalid_rows = 0;
    let mut matrix = if transpose { Matrix::new(col_count, row_count) } else { Matrix::new(row_count, col_count) };

    for (row_index, line) in lines[1..].iter().enumerate() {
        let items = hmf_common::java_text::split_all(line, ',');
        if items.len() != col_count + ignore_cols.len() {
            invalid_rows += 1;
            continue;
        }
        let mut col_index = 0;
        for (i, item) in items.iter().enumerate() {
            if ignore_cols.contains(&i) {
                continue;
            }
            let value = if item.is_empty() {
                0.0
            } else {
                hmf_common::java_text::parse_double(item)
                    .ok_or_else(|| anyhow!("java.lang.NumberFormatException: For input string: \"{item}\""))?
            };
            if transpose {
                matrix.data[col_index][row_index] = value;
            } else {
                matrix.data[row_index][col_index] = value;
            }
            col_index += 1;
        }
    }

    if invalid_rows > 0 {
        hmf_common::warn!("matrix(rows={} cols={}) has {} invalid rows", matrix.rows, matrix.cols, invalid_rows);
    }
    Ok(matrix)
}

/// MatrixFile.loadMatrixDataFile(filename, ...): the lines of the file (gzip if it ends in .gz), then as above;
/// Err where Java returns null (the file cannot be read, or holds no data rows)
pub fn load_matrix_file(
    filename: &str,
    column_names: &mut Vec<String>,
    ignore_fields: Option<&[&str]>,
    transpose: bool,
) -> Result<Matrix> {
    let lines = match hmf_common::java_text::read_lines(filename) {
        Ok(l) => l,
        Err(e) => {
            hmf_common::error!("failed to read matrix data file({filename}): {e}");
            bail!("failed to read matrix data file({filename}): {e}");
        }
    };
    let matrix = load_matrix_lines(&lines, column_names, ignore_fields, transpose)?;
    hmf_common::info!("loaded matrix(rows={} cols={}) from file({})", matrix.rows, matrix.cols, filename);
    Ok(matrix)
}
