//! End-to-end tests of the prep on a small BAM written here (synthetic reads): the comparisons with Java
//! are described in REPORT.md (synthetic inputs from scripts/make_synthetic.sh).

use rust_htslib::bam::{self, Read as BamRead};

const TELO: &str = "TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAG";
const PLAIN: &str = "ACGTACCTGATCGATCGTAGCTAGCTTGACGATCGATGCATGCAGTCAGTA";

struct Rec<'a> {
    name: &'a str,
    flag: u16,
    tid: i32,
    pos: i64,
    mtid: i32,
    mpos: i64,
    seq: &'a str,
    tags: Vec<(&'a [u8], String)>,
}

fn write_bam(path: &str, recs: &mut [Rec]) {
    let mut h = bam::Header::new();
    let mut hd = bam::header::HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate");
    h.push_record(&hd);
    for (n, l) in [("chr1", 1_000_000), ("chr2", 1_000_000)] {
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", n).push_tag(b"LN", l);
        h.push_record(&sq);
    }
    recs.sort_by_key(|r| if r.tid < 0 { (i32::MAX, 0) } else { (r.tid, r.pos) });
    {
        let mut w = bam::Writer::from_path(path, &h, bam::Format::Bam).unwrap();
        for r in recs.iter() {
            let mut b = bam::Record::new();
            let unmapped = r.flag & 4 != 0;
            let cigar = if unmapped { None } else { Some(bam::record::CigarString(vec![bam::record::Cigar::Match(r.seq.len() as u32)])) };
            b.set(r.name.as_bytes(), cigar.as_ref(), r.seq.as_bytes(), &vec![30u8; r.seq.len()]);
            b.set_flags(r.flag);
            b.set_tid(r.tid);
            b.set_pos(r.pos);
            b.set_mtid(r.mtid);
            b.set_mpos(r.mpos);
            b.set_mapq(if unmapped { 0 } else { 60 });
            for (t, v) in &r.tags {
                b.push_aux(t, bam::record::Aux::String(v.as_str())).unwrap();
            }
            w.write(&b).unwrap();
        }
    }
    bam::index::build(path, None, bam::index::Type::Bai, 1).unwrap();
}

fn names_in(path: &str) -> Vec<(String, u16)> {
    let mut r = bam::Reader::from_path(path).unwrap();
    let mut out = Vec::new();
    for rec in r.records() {
        let rec = rec.unwrap();
        out.push((String::from_utf8_lossy(rec.qname()).to_string(), rec.flags()));
    }
    out
}

#[test]
fn prep_groups_and_missing_reads() {
    let dir = std::env::temp_dir().join(format!("teal_e2e_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let bam_path = dir.join("in.bam").to_string_lossy().to_string();
    let telbam = dir.join("out.telbam.bam").to_string_lossy().to_string();
    let mut recs = vec![
        // a plain pair: not written
        Rec { name: "plain", flag: 99, tid: 0, pos: 1000, mtid: 0, mpos: 1200, seq: PLAIN, tags: vec![] },
        Rec { name: "plain", flag: 147, tid: 0, pos: 1200, mtid: 0, mpos: 1000, seq: PLAIN, tags: vec![] },
        // the mate comes first: found in the missing-read pass
        Rec { name: "late", flag: 99, tid: 0, pos: 5000, mtid: 1, mpos: 7000, seq: PLAIN, tags: vec![] },
        Rec { name: "late", flag: 147, tid: 1, pos: 7000, mtid: 0, mpos: 5000, seq: TELO, tags: vec![] },
        // the telomeric read first, its mate later and a supplementary read listed in SA
        Rec { name: "split", flag: 65, tid: 0, pos: 9000, mtid: 1, mpos: 9500, seq: TELO, tags: vec![(b"SA", format!("chr2,20001,+,{}M,60,0;", TELO.len()))] },
        Rec { name: "split", flag: 129, tid: 1, pos: 9500, mtid: 0, mpos: 9000, seq: PLAIN, tags: vec![] },
        Rec { name: "split", flag: 2113, tid: 1, pos: 20000, mtid: 1, mpos: 9500, seq: TELO, tags: vec![] },
        // a consensus read is skipped (and its mate with it)
        Rec { name: "cons", flag: 99, tid: 0, pos: 30000, mtid: 0, mpos: 30200, seq: TELO, tags: vec![(b"CR", "2;2;0".to_string())] },
        Rec { name: "cons", flag: 147, tid: 0, pos: 30200, mtid: 0, mpos: 30000, seq: PLAIN, tags: vec![] },
        // unplaced pair, the telomeric read second: its mate comes from the second look at the unplaced reads
        Rec { name: "unplaced", flag: 77, tid: -1, pos: -1, mtid: -1, mpos: -1, seq: PLAIN, tags: vec![] },
        Rec { name: "unplaced", flag: 141, tid: -1, pos: -1, mtid: -1, mpos: -1, seq: TELO, tags: vec![] },
        // a lone telomeric read whose mate is missing: written as incomplete
        Rec { name: "lone", flag: 97, tid: 1, pos: 40000, mtid: 1, mpos: 900000, seq: TELO, tags: vec![] },
        // a secondary alignment of a telomeric read: taken as supplementary
        Rec { name: "sec", flag: 65, tid: 0, pos: 50000, mtid: 0, mpos: 50300, seq: TELO, tags: vec![] },
        Rec { name: "sec", flag: 129, tid: 0, pos: 50300, mtid: 0, mpos: 50000, seq: PLAIN, tags: vec![] },
        // after the pair is complete: a new group, never complete
        Rec { name: "sec", flag: 321, tid: 1, pos: 60000, mtid: 0, mpos: 50300, seq: TELO, tags: vec![] },
    ];
    write_bam(&bam_path, &mut recs);
    for threads in [1, 2, 3] {
        crate::prep::process_bam(&crate::prep::TelbamParams { bam_file: &bam_path, ref_genome: None, telbam_file: &telbam, threads }).unwrap();
        let got = names_in(&telbam);
        let expected: Vec<(String, u16)> = [
            ("late", 99),
            ("split", 65),
            ("sec", 65),
            ("sec", 129),
            ("late", 147),
            ("split", 129),
            ("split", 2113),
            ("lone", 97),
            ("sec", 2113),
            ("unplaced", 77),
            ("unplaced", 141),
        ]
        .iter()
        .map(|(n, f)| (n.to_string(), *f))
        .collect();
        assert_eq!(got, expected, "threads {threads}");
    }
    let _ = std::fs::remove_dir_all(&dir);
}
