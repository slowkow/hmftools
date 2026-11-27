//! ViralInsertionTypeTest

use crate::prep::drivers::viral_insertion_type;

#[test]
fn test_from_virus_name() {
    for (name, t) in [
        ("Human papillomavirus type 16", "HPV"),
        ("Hepatitis B virus", "HBV"),
        ("Merkel cell polyomavirus", "MERKEL"),
        ("Human gammaherpesvirus 4", "HERPES"),
        ("Human gammaherpesvirus 8", "HERPES"),
    ] {
        assert_eq!(viral_insertion_type(name), t, "{name}");
    }
    assert_eq!(viral_insertion_type("Human mastadenovirus C"), "AAV");
    assert_eq!(viral_insertion_type("something else"), "OTHER");
}
