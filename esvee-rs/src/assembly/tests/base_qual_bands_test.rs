//! Not in Java: -low_base_qual_threshold and -high_base_qual_threshold (common::base_qual) move the bands
//! that assembly reads; the defaults are Java's constants.

use super::test_env;
use crate::assembly::assembly_utils::create_base_quals_above_min_threshold;
use crate::assembly::common_utils::{above_min_qual, below_min_qual, is_high_base_qual};
use crate::common::base_qual::set_test_base_qual_thresholds;

#[test]
fn base_qual_bands_follow_the_options() {
    let _env = test_env(); // Illumina, with the process-wide settings held
    assert!(below_min_qual(26) && above_min_qual(27));
    assert!(!is_high_base_qual(29) && is_high_base_qual(30));
    assert_eq!(create_base_quals_above_min_threshold(2), [27, 27]);

    set_test_base_qual_thresholds(Some((20, 24)));
    let moved = (below_min_qual(20), above_min_qual(21), is_high_base_qual(24), is_high_base_qual(23), create_base_quals_above_min_threshold(2));
    set_test_base_qual_thresholds(None);
    assert_eq!(moved, (true, true, true, false, vec![21, 21]));
}
