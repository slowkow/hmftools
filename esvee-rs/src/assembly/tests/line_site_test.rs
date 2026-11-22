//! LineSiteTest
use super::*;
use crate::assembly::common_utils::within_line_proximity;
use crate::assembly::constants::{LINE_BASE_A, LINE_BASE_T, LINE_INDEL_MAX_GAP, LINE_INDEL_MAX_OVERLAP};
use crate::assembly::junction_assembler::JunctionAssembler;
use crate::assembly::line_utils::{find_line_extension_end_index, try_line_sequence_link};
use crate::assembly::read::read_adjustments::mark_line_soft_clips;
use crate::assembly::types::StructuralVariantType;

#[test]
fn test_line_insertion_proximity() {
    let _env = test_env();
    assert!(!within_line_proximity(100, 100, FORWARD, FORWARD));
    assert!(within_line_proximity(100, 100, FORWARD, REVERSE));
    assert!(within_line_proximity(100, 100, REVERSE, FORWARD));
    // DEL within range
    let lower_position = 100;
    let upper_position = lower_position + LINE_INDEL_MAX_GAP - 1;
    assert!(within_line_proximity(lower_position, upper_position, FORWARD, REVERSE));
    assert!(within_line_proximity(upper_position, lower_position, REVERSE, FORWARD));
    // one further, still asserted as within range by the Java test
    assert!(within_line_proximity(lower_position, upper_position + 1, FORWARD, REVERSE));
    assert!(within_line_proximity(upper_position + 1, lower_position, REVERSE, FORWARD));
    // DUP within range
    let upper_position = lower_position + LINE_INDEL_MAX_OVERLAP - 1;
    assert!(within_line_proximity(lower_position, upper_position, REVERSE, FORWARD));
    assert!(within_line_proximity(upper_position, lower_position, FORWARD, REVERSE));
    assert!(within_line_proximity(lower_position, upper_position + 1, REVERSE, FORWARD));
    assert!(within_line_proximity(upper_position + 1, lower_position, FORWARD, REVERSE));
}

#[test]
fn test_read_line_extension_index() {
    let _env = test_env();
    let ref_bases = &REF_BASES_200[0..50];
    let line_sequence = "AAAAAAAAAAAAAAAA";
    let ext_bases = line_sequence.to_string();
    let read_bases = format!("{ext_bases}{ref_bases}");
    let read = create_read(&READ_ID_GENERATOR.next_id(), 100, &read_bases, &make_cigar_string(&read_bases, ext_bases.len(), 0));
    assert_eq!(0, find_line_extension_end_index(&read, LINE_BASE_A, ext_bases.len() as i32 - 1, false));

    // terminating after the min required length
    let ext_bases = format!("CCAAA{line_sequence}");
    let read_bases = format!("{ext_bases}{ref_bases}");
    let read = create_read(&READ_ID_GENERATOR.next_id(), 100, &read_bases, &make_cigar_string(&read_bases, ext_bases.len(), 0));
    assert_eq!(2, find_line_extension_end_index(&read, LINE_BASE_A, ext_bases.len() as i32 - 1, false));

    // forward junctions
    let line_sequence = reverse_complement_str("AAAAAAAAAAAAAAAA");
    let ext_bases = line_sequence.clone();
    let read_bases = format!("{ref_bases}{ext_bases}");
    let read = create_read(&READ_ID_GENERATOR.next_id(), 100, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len()));
    assert_eq!(read_bases.len() as i32 - 1, find_line_extension_end_index(&read, LINE_BASE_T, ref_bases.len() as i32, true));

    let ext_bases = format!("{line_sequence}TTTAA");
    let read_bases = format!("{ref_bases}{ext_bases}");
    let read = create_read(&READ_ID_GENERATOR.next_id(), 100, &read_bases, &make_cigar_string(&read_bases, 0, ext_bases.len()));
    assert_eq!(read_bases.len() as i32 - 3, find_line_extension_end_index(&read, LINE_BASE_T, ref_bases.len() as i32, true));
}

#[test]
fn test_line_junction_assembly() {
    let _env = test_env();
    let poly_at = "TGTTTATTTTTTTTTTTT"; // allows for 2 non-bases
    let extra_bases = "GTAGTGCTGT";
    let junction = new_junction(CHR_1, 100, FORWARD);

    // reads with poly-A tails of various lengths but all supporting the junction
    let ref_bases = &REF_BASES_400[51..101];
    let ext_bases1 = format!("{poly_at}{extra_bases}");
    let read_bases1 = format!("{ref_bases}{ext_bases1}");
    let read1 = create_read(&READ_ID_GENERATOR.next_id(), 51, &read_bases1, &make_cigar_string(&read_bases1, 0, ext_bases1.len()));
    let read1b = clone_read(&read1, &READ_ID_GENERATOR.next_id());
    let ext_bases2 = format!("{poly_at}TT{extra_bases}");
    let read_bases2 = format!("{}{}", &ref_bases[1..], ext_bases2);
    let read2 = create_read(&READ_ID_GENERATOR.next_id(), 52, &read_bases2, &make_cigar_string(&read_bases2, 0, ext_bases2.len()));
    let ext_bases3 = format!("{poly_at}T{extra_bases}");
    let read_bases3 = format!("{ref_bases}{ext_bases3}");
    let read3 = create_read(&READ_ID_GENERATOR.next_id(), 51, &read_bases3, &make_cigar_string(&read_bases3, 0, ext_bases3.len()));
    // too short to be marked as line but supports the junction
    let ext_bases4 = &poly_at[0..12];
    let read_bases4 = format!("{ref_bases}{ext_bases4}");
    let read4 = create_read(&READ_ID_GENERATOR.next_id(), 51, &read_bases4, &make_cigar_string(&read_bases4, 0, ext_bases4.len()));

    let reads = vec![read1.clone(), read1b, read2.clone(), read3.clone(), read4.clone()];
    for read in &reads {
        mark_line_soft_clips(read.m());
    }
    assert!(read1.has_line_tail_end());
    assert!(read2.has_line_tail_end());
    assert!(read3.has_line_tail_end());
    assert!(!read1.has_line_tail_start());
    assert!(!read4.has_line_tail_start());
    assert!(!read4.has_line_tail_end());

    let mut junction_assembler = JunctionAssembler::new(junction, None, None);
    let assemblies = junction_assembler.process_junction(&reads);
    assert_eq!(1, assemblies.len());
    let assembly = assemblies[0];
    assert!(assembly.has_line_sequence());
    assert_eq!(5, assembly.support_count());
    assert_eq!(ext_bases1, assembly.form_junction_sequence()); // based on taking the median T length
}

#[test]
fn test_line_assembly_links() {
    let _env = test_env();
    let first_ref_bases = &REF_BASES_200[0..100];
    let second_ref_bases = &REF_BASES_200[110..200];
    let poly_a = "AAAAAAAAAAAAAAAA";
    let poly_t = reverse_complement_str(poly_a);

    // first a DEL with poly A at -ve junction
    let first_junction = new_junction(CHR_1, 100, FORWARD);
    let second_junction = new_junction(CHR_1, 110, REVERSE);
    let extra_bases = "GTAGTGCTGTCGC";
    let first_ext_bases = format!("{extra_bases}{poly_a}");
    let first_assembly_bases = format!("{first_ref_bases}{first_ext_bases}");
    let base_quals = build_default_base_quals(first_assembly_bases.len());
    let first_assembly = JunctionAssembly::new_test(first_junction, first_assembly_bases.as_bytes(), &base_quals, first_ref_bases.len() as i32 - 1);
    let second_ext_bases = poly_a.to_string();
    let second_assembly_bases = format!("{second_ext_bases}{second_ref_bases}");
    let second_assembly = JunctionAssembly::new_test(second_junction, second_assembly_bases.as_bytes(), &base_quals, second_ext_bases.len() as i32);
    second_assembly.m().mark_line_sequence();
    let link = try_line_sequence_link(first_assembly, second_assembly, false, false).expect("link");
    assert_eq!(first_ext_bases, link.inserted_bases());
    assert!(first_assembly == link.first());
    assert!(second_assembly == link.second());
    assert_eq!(StructuralVariantType::Del, link.sv_type());

    // then a DEL with poly T at +ve junction
    let first_ext_bases = poly_t.clone();
    let first_assembly_bases = format!("{first_ref_bases}{first_ext_bases}");
    let base_quals = build_default_base_quals(first_assembly_bases.len());
    let first_assembly = JunctionAssembly::new_test(first_junction, first_assembly_bases.as_bytes(), &base_quals, first_ref_bases.len() as i32 - 1);
    first_assembly.m().mark_line_sequence();
    let second_ext_bases = format!("{poly_t}{extra_bases}");
    let second_assembly_bases = format!("{second_ext_bases}{second_ref_bases}");
    let second_assembly = JunctionAssembly::new_test(second_junction, second_assembly_bases.as_bytes(), &base_quals, second_ext_bases.len() as i32);
    let link = try_line_sequence_link(first_assembly, second_assembly, false, false).expect("link");
    assert_eq!(second_ext_bases, link.inserted_bases());
    assert!(first_assembly == link.first());
    assert!(second_assembly == link.second());
    assert_eq!(StructuralVariantType::Del, link.sv_type());

    // same orientation
    let first_ext_bases = format!("{poly_t}{extra_bases}");
    let first_assembly_bases = format!("{first_ext_bases}{first_ref_bases}");
    let first_junction = new_junction(CHR_1, 100, REVERSE);
    let first_assembly = JunctionAssembly::new_test(first_junction, first_assembly_bases.as_bytes(), &base_quals, first_ext_bases.len() as i32);
    let second_ext_bases = poly_a.to_string();
    let second_assembly_bases = format!("{second_ext_bases}{second_ref_bases}");
    let second_assembly = JunctionAssembly::new_test(second_junction, second_assembly_bases.as_bytes(), &base_quals, second_ext_bases.len() as i32);
    second_assembly.m().mark_line_sequence();
    let link = try_line_sequence_link(first_assembly, second_assembly, true, false).expect("link");
    assert_eq!(reverse_complement_str(&first_ext_bases), link.inserted_bases());
    assert_eq!(StructuralVariantType::Inv, link.sv_type());

    // first assembly has the poly-T sequence
    let extra_bases2 = "GTCGT";
    let first_ext_bases = format!("{poly_t}{extra_bases2}{poly_t}");
    let first_assembly_bases = format!("{first_ref_bases}{first_ext_bases}");
    let first_junction = new_junction(CHR_1, 120, FORWARD);
    let first_assembly = JunctionAssembly::new_test(first_junction, first_assembly_bases.as_bytes(), &base_quals, first_ref_bases.len() as i32 - 1);
    first_assembly.m().mark_line_sequence();
    let second_ext_bases = format!("{poly_t}{extra_bases2}{poly_t}{extra_bases}");
    let second_assembly_bases = format!("{second_ext_bases}{second_ref_bases}");
    let second_assembly = JunctionAssembly::new_test(second_junction, second_assembly_bases.as_bytes(), &base_quals, second_ext_bases.len() as i32);
    let link = try_line_sequence_link(first_assembly, second_assembly, false, false).expect("link");
    assert_eq!(format!("{first_ext_bases}{extra_bases}"), link.inserted_bases());
    assert_eq!(StructuralVariantType::Dup, link.sv_type());
}
