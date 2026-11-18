//! Cursor-based selectors: PanelSelector, HotspotCache, TierSelector. The cursors are stateful and
//! path-dependent in Java, so they are ported literally.

use crate::types::*;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PanelStatus {
    WithinPanel,
    OutsidePanel,
    Mixed,
}

pub struct PanelSelector<'a> {
    regions: &'a [BaseRegion],
    index: usize,
}

impl<'a> PanelSelector<'a> {
    pub fn new(regions: &'a [BaseRegion]) -> Self {
        PanelSelector { regions, index: 0 }
    }
    fn update(&mut self, start: i32, end: i32) {
        if self.regions.is_empty() {
            return;
        }
        while self.index > 0 && self.regions[self.index].start > end {
            self.index -= 1;
        }
        while self.index < self.regions.len() - 1 && self.regions[self.index].end < start {
            self.index += 1;
        }
    }
    pub fn panel_status(&mut self, start: i32, end: i32) -> PanelStatus {
        if self.regions.is_empty() {
            return PanelStatus::OutsidePanel;
        }
        self.update(start, end);
        let c = self.regions[self.index];
        if positions_within(start, end, c.start, c.end) {
            return PanelStatus::WithinPanel;
        }
        if start > c.end || end < c.start {
            return PanelStatus::OutsidePanel;
        }
        PanelStatus::Mixed
    }
    pub fn in_panel(&mut self, start: i32, end: i32) -> bool {
        if self.regions.is_empty() {
            return false;
        }
        self.update(start, end);
        let c = self.regions[self.index];
        start <= c.end && end >= c.start
    }
}

pub const NEAR_HOTSPOT_DISTANCE: i32 = 5;

pub struct HotspotCache<'a> {
    hotspots: &'a [SimpleVariant],
    min_pos: i32,
    max_pos: i32,
    index: i64,
}

impl<'a> HotspotCache<'a> {
    pub fn new(hotspots: &'a [SimpleVariant]) -> Self {
        let (min_pos, max_pos) = match (hotspots.first(), hotspots.last()) {
            (Some(f), Some(l)) => (f.pos, l.pos),
            _ => (-1, -1),
        };
        HotspotCache { hotspots, min_pos, max_pos, index: 0 }
    }
    /// HotspotCache.setIndex (@VisibleForTesting in Java)
    #[cfg(test)]
    pub(crate) fn set_index(&mut self, index: i64) {
        self.index = index.min(self.hotspots.len() as i64 - 1).max(0);
    }
    fn at(&self, i: i64) -> &SimpleVariant {
        &self.hotspots[i as usize]
    }
    pub fn is_hotspot(&mut self, v: &SimpleVariant) -> bool {
        if self.hotspots.is_empty() {
            return false;
        }
        if !position_within(v.pos, self.min_pos, self.max_pos) {
            return false;
        }
        let n = self.hotspots.len() as i64;
        while self.index >= 0 && self.at(self.index).pos == v.pos {
            self.index -= 1;
            if self.index < 0 || self.at(self.index).pos < v.pos {
                self.index += 1;
                break;
            }
        }
        if self.at(self.index).pos <= v.pos {
            if self.at(self.index) == v {
                return true;
            }
            self.index += 1;
            while self.index < n {
                let h = self.at(self.index);
                if h.pos > v.pos {
                    break;
                }
                if h == v {
                    return true;
                }
                self.index += 1;
            }
        } else {
            self.index -= 1;
            while self.index >= 0 {
                let h = self.at(self.index);
                if h.pos < v.pos {
                    return false;
                }
                if h == v {
                    return true;
                }
                self.index -= 1;
            }
        }
        self.index = self.index.max(0).min(n - 1);
        false
    }

    fn is_near(h: &SimpleVariant, start: i32, end: i32) -> bool {
        positions_overlap(h.pos - NEAR_HOTSPOT_DISTANCE, h.pos + h.ref_len() - 1 + NEAR_HOTSPOT_DISTANCE, start, end)
    }

    pub fn is_near_hotspot(&mut self, v: &SimpleVariant) -> bool {
        if self.hotspots.is_empty() {
            return false;
        }
        let n = self.hotspots.len() as i64;
        let start = v.pos;
        let end = v.position_end();
        if self.at(self.index).pos <= start {
            if Self::is_near(self.at(self.index), start, end) {
                return true;
            }
            self.index += 1;
            while self.index < n {
                let h = self.at(self.index);
                if Self::is_near(h, start, end) {
                    return true;
                }
                if h.pos > start {
                    return false;
                }
                self.index += 1;
            }
        } else {
            if Self::is_near(self.at(self.index), start, end) {
                return true;
            }
            self.index -= 1;
            while self.index >= 0 {
                let h = self.at(self.index);
                if Self::is_near(h, start, end) {
                    return true;
                }
                if h.pos < start {
                    return false;
                }
                self.index -= 1;
            }
        }
        self.index = self.index.max(0).min(n - 1);
        false
    }
}

pub struct TierSelector<'a> {
    hotspots: HotspotCache<'a>,
    panel: PanelSelector<'a>,
    high_conf: PanelSelector<'a>,
}

impl<'a> TierSelector<'a> {
    pub fn new(hotspots: &'a [SimpleVariant], panel: &'a [BaseRegion], high_conf: &'a [BaseRegion]) -> Self {
        TierSelector {
            hotspots: HotspotCache::new(hotspots),
            panel: PanelSelector::new(panel),
            high_conf: PanelSelector::new(high_conf),
        }
    }
    pub fn tier(&mut self, v: &SimpleVariant) -> Tier {
        if self.hotspots.is_hotspot(v) {
            Tier::Hotspot
        } else if self.panel.in_panel(v.pos, v.end()) {
            Tier::Panel
        } else if self.high_conf.in_panel(v.pos, v.end()) {
            Tier::HighConfidence
        } else {
            Tier::LowConfidence
        }
    }
}
