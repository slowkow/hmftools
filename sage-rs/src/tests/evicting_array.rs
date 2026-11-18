//! Ported from sage/src/test/java/com/hartwig/hmftools/sage/common/EvictingArrayTest.java

use super::util::*;
use crate::candidate::{build_region_blocks, EvictingArray, RefContext};
use crate::select::{PanelSelector, PanelStatus};
use crate::types::*;

// Java: testRegionBlocks
#[test]
fn test_region_blocks() {
    let panel = vec![BaseRegion::new(950, 1050), BaseRegion::new(1100, 1450)];

    let chromosome = "6";

    // HlaCommon.populateGeneData with one HLA gene at 1000-1300 on chromosome 6
    let mut ref_data = empty_ref_data();
    ref_data.hla_chromosome = chromosome.to_string();
    ref_data.hla_genes = vec![BaseRegion::new(1000, 1300)];

    let region = BaseRegion::new(700, 1500);

    let hotspots = vec![SimpleVariant::new(910, b"A", b"T"), SimpleVariant::new(1350, b"A", b"T")];

    // the block size is REGION_BLOCK_SIZE inside build_region_blocks
    let mut region_blocks = build_region_blocks(chromosome, region, &panel, &hotspots, &ref_data, 20, 10);

    assert_eq!(8, region_blocks.len());

    assert_eq!(10, region_blocks[0].depth_limit);
    assert_eq!(10, region_blocks[1].depth_limit);
    assert_eq!(20, region_blocks[2].depth_limit);
    // the Rust RegionBlock does not keep its panel status (only the depth limit derived from it), so
    // RegionBlock.panelStatus() is recomputed here the way build_region_blocks does
    let mut panel_selector = PanelSelector::new(&panel);
    assert_eq!(PanelStatus::Mixed, panel_selector.panel_status(region_blocks[2].start, region_blocks[2].end));
    assert_eq!(20, region_blocks[7].depth_limit);

    assert!(region_blocks[0].apply_event_penalty);
    assert!(!region_blocks[4].apply_event_penalty);

    assert!(region_blocks[2].covers_hotspot(900, 920));

    let region_block = &mut region_blocks[0];

    for _ in 0..9 {
        region_block.depth += 1; // incrementDepth
    }

    assert!(!region_block.limit_reached());
    region_block.depth += 1;
    assert!(region_block.limit_reached());
}

// Java: testEvictingArray
#[test]
fn test_evicting_array() {
    let read_length = 200;

    // the eviction handler is the vector the array appends evicted contexts to
    let mut handler_items: Vec<RefContext> = Vec::new();
    let mut evicting_array = EvictingArray::with_capacity(400, 280);

    let capacity = evicting_array.capacity();
    assert_eq!(read_length * 2, capacity);

    let start_position = 1000;
    let min_position = start_position - evicting_array.read_length_buffer();
    evicting_array.ensure(start_position, &mut handler_items);
    assert_eq!(min_position, evicting_array.min_position());

    let mut position = start_position;
    while position < min_position + capacity {
        evicting_array.ensure(position, &mut handler_items);
        position += 1;
    }

    assert_eq!(120, evicting_array.item_count());

    // no flush yet
    assert_eq!(0, handler_items.len());

    while position < 1400 {
        evicting_array.ensure(position, &mut handler_items);
        position += 1;
    }

    assert_eq!(1080, evicting_array.min_position());
    assert_eq!(320, evicting_array.item_count());
    assert_eq!(80, handler_items.len());

    while position < 1600 {
        evicting_array.ensure(position, &mut handler_items);
        position += 1;
    }

    assert_eq!(1200, evicting_array.min_position());
    assert_eq!(400, evicting_array.item_count());
    assert_eq!(200, handler_items.len());

    // evict all but the last

    evicting_array.ensure(2000, &mut handler_items);
    assert_eq!(1, evicting_array.item_count());
    assert_eq!(600, handler_items.len());

    // evictAll() is checkFlush(-1)
    evicting_array.check_flush(-1, &mut handler_items);
    assert_eq!(0, evicting_array.item_count());
    assert_eq!(601, handler_items.len());
}
