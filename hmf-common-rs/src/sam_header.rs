//! SAM header lines as htsjdk writes them.
//!
//! htsjdk's text header decoder does not keep the input's `@HD VN`: the decoded SAMFileHeader keeps the
//! version its constructor sets, SAMFileHeader.CURRENT_VERSION (1.6 in htsjdk 2.23), and the BAM writer
//! writes that. So a Java tool that writes a BAM with the header of the BAM it read (REDUX, TEAL's
//! telbam, CIDER's cider.bam) writes `VN:1.6` whatever the input says: `VN:1.3` from bwa-mem3 and `VN:1.5`
//! from other aligners become `VN:1.6`. The other tags of the line are kept in their order.

/// SAMFileHeader.CURRENT_VERSION of htsjdk 2.23
pub const HTSJDK_SAM_VERSION: &str = "1.6";

/// The `@HD` line `line` as htsjdk writes it: `VN:1.6` first, then the line's other tags in their order
pub fn htsjdk_hd_line(line: &str) -> String {
    let mut out = format!("@HD\tVN:{HTSJDK_SAM_VERSION}");
    for f in line.split('\t').skip(1).filter(|f| !f.starts_with("VN:") && !f.is_empty()) {
        out.push('\t');
        out.push_str(f);
    }
    out
}

/// The header text `text` with its `@HD` line (if any) as htsjdk writes it; the other lines unchanged
pub fn htsjdk_header_version(text: &str) -> String {
    let mut out = String::with_capacity(text.len() + 4);
    for line in text.split_inclusive('\n') {
        let (body, nl) = match line.strip_suffix('\n') {
            Some(b) => (b, "\n"),
            None => (line, ""),
        };
        if body.starts_with("@HD\t") || body == "@HD" {
            out.push_str(&htsjdk_hd_line(body));
        } else {
            out.push_str(body);
        }
        out.push_str(nl);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hd_line_takes_htsjdk_version() {
        assert_eq!(htsjdk_hd_line("@HD\tVN:1.3\tSO:coordinate"), "@HD\tVN:1.6\tSO:coordinate");
        assert_eq!(htsjdk_hd_line("@HD\tVN:1.6\tSO:coordinate"), "@HD\tVN:1.6\tSO:coordinate");
        assert_eq!(htsjdk_hd_line("@HD\tSO:queryname\tVN:1.5\tGO:query"), "@HD\tVN:1.6\tSO:queryname\tGO:query");
        assert_eq!(htsjdk_hd_line("@HD"), "@HD\tVN:1.6");
    }

    #[test]
    fn header_text_changes_only_hd() {
        let t = "@HD\tVN:1.5\tSO:coordinate\n@SQ\tSN:chr1\tLN:10\n@CO\tVN:1.5 stays\n";
        assert_eq!(htsjdk_header_version(t), "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10\n@CO\tVN:1.5 stays\n");
        assert_eq!(htsjdk_header_version("@SQ\tSN:a\tLN:1"), "@SQ\tSN:a\tLN:1");
    }
}
