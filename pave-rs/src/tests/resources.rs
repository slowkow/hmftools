//! ResourcesTest

use super::*;
use hmf_common::vcf::VcfRecord;

use crate::annotation::clinvar::ClinvarData;
use crate::constants::CLNSIG;

/// createVariant: a variant with an empty record
pub fn create_variant(chromosome: &str, position: i32, r: &str, a: &str) -> VariantData {
    let mut v = VariantData::new(chromosome, position, r, a);
    v.set_context(VcfRecord::new(chromosome, position, ".", r, a));
    v
}

#[test]
fn test_clinvar() {
    let mut var1 = create_variant(CHR_1, 100, "A", "G");
    let mut var2 = create_variant(CHR_1, 100, "A", "C");
    let mut var3 = create_variant(CHR_1, 100, "A", "T");

    let mut cache = ClinvarData::default();
    let (sig, conflict) = ("pathogenic", "conflict");
    cache.add_entry(CHR_1, 99, "A", "T", sig, conflict);
    cache.add_entry(CHR_1, 100, "A", "AT", sig, conflict);
    cache.add_entry(CHR_1, 100, "A", "C", sig, conflict);
    cache.add_entry(CHR_1, 100, "A", "G", sig, conflict);
    cache.add_entry(CHR_1, 101, "A", "T", sig, conflict);

    for v in [&mut var1, &mut var2, &mut var3] {
        let (pos, r, a) = (v.position, v.ref_.clone(), v.alt.clone());
        cache.annotate(v.context_mut(), CHR_1, pos, &r, &a).unwrap();
    }
    assert!(var1.context().has_attr(CLNSIG));
    assert!(var2.context().has_attr(CLNSIG));
    assert!(!var3.context().has_attr(CLNSIG));
}
