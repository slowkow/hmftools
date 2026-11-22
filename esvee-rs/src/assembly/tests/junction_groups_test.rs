//! JunctionGroupsTest
use std::collections::HashMap;

use super::*;
use crate::assembly::types::junction::JuncP;

#[test]
fn test_junction_group_merging() {
    let _env = test_env();
    let mut existing_junctions: Vec<JuncP> = Vec::new();
    let mut new_junctions: Vec<JuncP> = Vec::new();
    for (pos, orient) in [(100, FORWARD), (110, FORWARD), (120, FORWARD), (120, REVERSE), (130, REVERSE), (140, REVERSE)] {
        existing_junctions.push(new_junction(CHR_1, pos, orient));
    }
    // 120 forward, 120 reverse and 140 reverse are duplicates
    for (pos, orient) in [(95, FORWARD), (115, FORWARD), (120, FORWARD), (120, REVERSE), (130, FORWARD), (135, REVERSE), (140, REVERSE), (145, REVERSE)] {
        new_junctions.push(new_junction(CHR_1, pos, orient));
    }
    let mut existing_map: HashMap<String, Vec<JuncP>> = HashMap::new();
    existing_map.insert(CHR_1.to_string(), existing_junctions);
    let mut new_map: HashMap<String, Vec<JuncP>> = HashMap::new();
    new_map.insert(CHR_1.to_string(), new_junctions);

    Junction::merge_junctions(&mut existing_map, &new_map);

    assert_eq!(11, existing_map[CHR_1].len());
    assert!(Junction::validate_junction_map(&existing_map));
}
