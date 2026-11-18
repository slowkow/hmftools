//! Ported from sage/src/test/java/com/hartwig/hmftools/sage/common/VariantSelectorsTest.java
//!
//! The Rust SimpleVariant has no chromosome (it is implicit), so `new SimpleVariant(CHR_1, pos, ref, alt)`
//! becomes `SimpleVariant::new(pos, ref, alt)`.

use crate::select::{HotspotCache, PanelSelector, PanelStatus};
use crate::types::*;

fn panel() -> Vec<BaseRegion> {
    vec![BaseRegion::new(995, 995), BaseRegion::new(998, 1102), BaseRegion::new(1995, 1995), BaseRegion::new(1998, 2102)]
}

fn check_overlap(panel_selector: &mut PanelSelector) {
    assert!(panel_selector.in_panel(995, 995));
    assert!(panel_selector.in_panel(994, 996));
    assert!(panel_selector.in_panel(1102, 1102));
    assert!(panel_selector.in_panel(1000, 1002));
    assert!(panel_selector.in_panel(996, 1000));
    assert!(!panel_selector.in_panel(1, 994));
    assert!(panel_selector.in_panel(998, 998));
    assert!(!panel_selector.in_panel(996, 997));
    assert!(panel_selector.in_panel(994, 1000));
    assert!(!panel_selector.in_panel(2200, 10000));
}

// Java: testOverlap
#[test]
fn test_overlap() {
    let panel = panel();
    let mut panel_selector = PanelSelector::new(&panel);
    check_overlap(&mut panel_selector);
}

// Java: testOutOfOrder
#[test]
fn test_out_of_order() {
    let panel = panel();
    let mut panel_selector = PanelSelector::new(&panel);
    check_overlap(&mut panel_selector);
    check_overlap(&mut panel_selector);
}

// Java: testPanelStatus
#[test]
fn test_panel_status() {
    let panel = panel();
    let mut panel_selector = PanelSelector::new(&panel);

    // panelStatus(position) is panelStatus(position, position)
    assert_eq!(PanelStatus::WithinPanel, panel_selector.panel_status(995, 995));
    assert_eq!(PanelStatus::OutsidePanel, panel_selector.panel_status(996, 996));
    assert_eq!(PanelStatus::OutsidePanel, panel_selector.panel_status(994, 994));
    assert_eq!(PanelStatus::OutsidePanel, panel_selector.panel_status(2103, 2103));
    assert_eq!(PanelStatus::WithinPanel, panel_selector.panel_status(2100, 2100));

    assert_eq!(PanelStatus::Mixed, panel_selector.panel_status(995, 996));
    assert_eq!(PanelStatus::Mixed, panel_selector.panel_status(994, 995));
    assert_eq!(PanelStatus::Mixed, panel_selector.panel_status(994, 996)); // straddling
    assert_eq!(PanelStatus::Mixed, panel_selector.panel_status(1995, 2000));
    assert_eq!(PanelStatus::Mixed, panel_selector.panel_status(1500, 2500)); // straddling
}

// Java: testHotspotSelection
#[test]
fn test_hotspot_selection() {
    let r = b"A";
    let hotspots = vec![
        SimpleVariant::new(100, r, b"C"),
        SimpleVariant::new(101, r, b"C"),
        SimpleVariant::new(101, r, b"G"),
        SimpleVariant::new(101, r, b"T"),
        SimpleVariant::new(102, r, b"C"),
        SimpleVariant::new(103, r, b"C"),
        SimpleVariant::new(103, r, b"G"),
        SimpleVariant::new(103, r, b"T"),
        SimpleVariant::new(104, r, b"C"),
    ];

    let mut hotspot_cache = HotspotCache::new(&hotspots);

    let mut test_variant = SimpleVariant::new(99, r, b"C");
    assert!(!hotspot_cache.is_hotspot(&test_variant));

    test_variant = SimpleVariant::new(105, r, b"C");
    assert!(!hotspot_cache.is_hotspot(&test_variant));

    test_variant = SimpleVariant::new(100, r, b"G");
    assert!(!hotspot_cache.is_hotspot(&test_variant));

    test_variant = SimpleVariant::new(100, r, b"C");
    assert!(hotspot_cache.is_hotspot(&test_variant));

    test_variant = SimpleVariant::new(101, r, b"C");
    assert!(hotspot_cache.is_hotspot(&test_variant));

    test_variant = SimpleVariant::new(101, r, b"T");
    assert!(hotspot_cache.is_hotspot(&test_variant));

    // moving backwards through matches
    test_variant = SimpleVariant::new(101, r, b"G");
    assert!(hotspot_cache.is_hotspot(&test_variant));

    test_variant = SimpleVariant::new(104, r, b"G");
    assert!(!hotspot_cache.is_hotspot(&test_variant));

    test_variant = SimpleVariant::new(104, r, b"C");
    assert!(hotspot_cache.is_hotspot(&test_variant));

    test_variant = SimpleVariant::new(103, r, b"C");
    assert!(hotspot_cache.is_hotspot(&test_variant));

    test_variant = SimpleVariant::new(103, r, b"T");
    assert!(hotspot_cache.is_hotspot(&test_variant));

    test_variant = SimpleVariant::new(103, r, b"G");
    assert!(hotspot_cache.is_hotspot(&test_variant));

    test_variant = SimpleVariant::new(102, r, b"G");
    assert!(!hotspot_cache.is_hotspot(&test_variant));

    test_variant = SimpleVariant::new(102, r, b"C");
    assert!(hotspot_cache.is_hotspot(&test_variant));
}

// Java: testHotspotProximity
#[test]
fn test_hotspot_proximity() {
    let r = b"A";
    let alt = b"C";

    let hotspots = vec![
        SimpleVariant::new(100, r, alt),
        SimpleVariant::new(110, r, alt),
        SimpleVariant::new(110, r, alt),
        SimpleVariant::new(122, r, alt),
        SimpleVariant::new(123, r, alt),
        SimpleVariant::new(123, r, alt),
        SimpleVariant::new(123, b"AAA", alt),
        SimpleVariant::new(140, r, alt),
    ];

    let mut hotspot_cache = HotspotCache::new(&hotspots);

    let mut test_variant = SimpleVariant::new(94, r, alt);
    assert!(!hotspot_cache.is_near_hotspot(&test_variant));

    test_variant = SimpleVariant::new(95, r, alt);
    assert!(hotspot_cache.is_near_hotspot(&test_variant));

    test_variant = SimpleVariant::new(105, r, alt);
    assert!(hotspot_cache.is_near_hotspot(&test_variant));

    test_variant = SimpleVariant::new(116, r, alt);
    assert!(!hotspot_cache.is_near_hotspot(&test_variant));

    // test navigation
    test_variant = SimpleVariant::new(124, r, alt);
    assert!(hotspot_cache.is_near_hotspot(&test_variant));

    test_variant = SimpleVariant::new(120, r, alt);
    assert!(hotspot_cache.is_near_hotspot(&test_variant));

    test_variant = SimpleVariant::new(135, r, alt);
    assert!(hotspot_cache.is_near_hotspot(&test_variant));

    test_variant = SimpleVariant::new(145, r, alt);
    assert!(hotspot_cache.is_near_hotspot(&test_variant));

    test_variant = SimpleVariant::new(146, r, alt);
    assert!(!hotspot_cache.is_near_hotspot(&test_variant));

    hotspot_cache.set_index(0);

    test_variant = SimpleVariant::new(124, r, alt);
    assert!(hotspot_cache.is_near_hotspot(&test_variant));

    hotspot_cache.set_index(hotspots.len() as i64 - 1);

    test_variant = SimpleVariant::new(119, r, alt);
    assert!(hotspot_cache.is_near_hotspot(&test_variant));

    // check indels
    test_variant = SimpleVariant::new(130, r, alt);
    assert!(hotspot_cache.is_near_hotspot(&test_variant));

    test_variant = SimpleVariant::new(133, b"AAA", alt);
    assert!(hotspot_cache.is_near_hotspot(&test_variant));
}
