//! Memory diagnostics (not part of the Java code): counts of live objects of the kinds that hold most of the
//! memory, and the process memory footprint at stage boundaries. The counters are always kept (one relaxed
//! atomic operation per object); the log lines are only written when the environment variable
//! `ESVEE_MEM_STATS` is set, so normal log output is unchanged.

use std::sync::atomic::{AtomicBool, AtomicI64, Ordering::Relaxed};

pub static READS_CREATED: AtomicI64 = AtomicI64::new(0);
pub static READS_LIVE: AtomicI64 = AtomicI64::new(0);
/// estimated heap bytes of the live reads (struct plus the buffers it owns at creation)
pub static READ_BYTES_LIVE: AtomicI64 = AtomicI64::new(0);
pub static SUPPORT_READS_LIVE: AtomicI64 = AtomicI64::new(0);
pub static ASSEMBLIES_CREATED: AtomicI64 = AtomicI64::new(0);
pub static ASSEMBLIES_RELEASED: AtomicI64 = AtomicI64::new(0);
/// reads that were still referenced from elsewhere when their junction group was done
pub static READS_RETAINED: AtomicI64 = AtomicI64::new(0);

static ENABLED: AtomicBool = AtomicBool::new(false);

pub fn init_from_env() {
    ENABLED.store(std::env::var_os("ESVEE_MEM_STATS").is_some(), Relaxed);
}
pub fn enabled() -> bool {
    ENABLED.load(Relaxed)
}

/// Current memory footprint of the process in bytes (the figure `/usr/bin/time -l` reports the peak of),
/// or zero where it is not available.
#[cfg(target_os = "macos")]
pub fn memory_footprint() -> u64 {
    unsafe extern "C" {
        fn proc_pid_rusage(pid: i32, flavor: i32, buffer: *mut u64) -> i32;
    }
    // rusage_info_v0: a 16 byte uuid followed by 64-bit counters, the eighth being the footprint
    let mut buffer = [0u64; 12];
    let rc = unsafe { proc_pid_rusage(std::process::id() as i32, 0, buffer.as_mut_ptr()) };
    if rc == 0 { buffer[9] } else { 0 }
}
#[cfg(not(target_os = "macos"))]
pub fn memory_footprint() -> u64 {
    0
}

pub fn log_stage(stage: &str) {
    if !enabled() {
        return;
    }
    const MB: f64 = 1024.0 * 1024.0;
    crate::info!(
        "memory at {}: footprint({:.0}MB) reads(live={} created={} retained={} approx={:.0}MB) supportReads({}) assemblies(created={} released={})",
        stage,
        memory_footprint() as f64 / MB,
        READS_LIVE.load(Relaxed),
        READS_CREATED.load(Relaxed),
        READS_RETAINED.load(Relaxed),
        READ_BYTES_LIVE.load(Relaxed) as f64 / MB,
        SUPPORT_READS_LIVE.load(Relaxed),
        ASSEMBLIES_CREATED.load(Relaxed),
        ASSEMBLIES_RELEASED.load(Relaxed),
    );
    log_allocator_size_classes(stage);
}

/// Logs an estimate of what the assemblies of the junction groups hold (support reads counted once).
pub fn log_assembly_breakdown(junction_groups: &[crate::assembly::types::junction_group::JuncGroupP]) {
    if !enabled() {
        return;
    }
    use std::collections::HashSet;
    let mut seen_support: HashSet<usize> = HashSet::new();
    let mut seen_regions: HashSet<usize> = HashSet::new();
    let (mut assemblies, mut assembly_bytes, mut support, mut support_bytes) = (0usize, 0usize, 0usize, 0usize);
    let (mut regions, mut region_ids, mut region_id_bytes, mut soft_clip_ids) = (0usize, 0usize, 0usize, 0usize);
    for junction_group in junction_groups {
        for assembly in junction_group.junction_assemblies() {
            assemblies += 1;
            assembly_bytes += std::mem::size_of::<crate::assembly::types::junction_assembly::JunctionAssembly>()
                + assembly.bases().len()
                + assembly.base_quals().len()
                + assembly.support().capacity() * 8;
            for s in assembly.support() {
                if seen_support.insert(s.addr()) {
                    support += 1;
                    support_bytes += s.approx_bytes();
                }
            }
            for region in assembly.remote_regions() {
                if seen_regions.insert(region.addr()) {
                    regions += 1;
                    for id in region.read_ids().iter_unordered() {
                        region_ids += 1;
                        region_id_bytes += id.capacity() + 24;
                    }
                }
            }
            for soft_clip in assembly.ref_side_soft_clips() {
                soft_clip_ids += soft_clip.read_ids().len();
            }
        }
    }
    const MB: f64 = 1024.0 * 1024.0;
    crate::info!(
        "memory held by junction assemblies: assemblies({} approx={:.0}MB) supportReads({} approx={:.0}MB structSize={}) remoteRegions({} readIds={} approx={:.0}MB) refSideSoftClipReadIds({}) sizes(read={} assembly={})",
        assemblies,
        assembly_bytes as f64 / MB,
        support,
        support_bytes as f64 / MB,
        std::mem::size_of::<crate::assembly::types::support_read::SupportRead>(),
        regions,
        region_ids,
        region_id_bytes as f64 / MB,
        soft_clip_ids,
        std::mem::size_of::<crate::assembly::read::read::Read>(),
        std::mem::size_of::<crate::assembly::types::junction_assembly::JunctionAssembly>(),
    );
}

/// Per block size of the allocator: committed bytes and bytes in use (only with the mimalloc allocator).
#[cfg(feature = "mimalloc")]
pub fn log_allocator_size_classes(stage: &str) {
    if !enabled() {
        return;
    }
    use std::collections::BTreeMap;
    use std::ffi::c_void;
    #[repr(C)]
    struct Area {
        blocks: *mut c_void,
        reserved: usize,
        committed: usize,
        used: usize,
        block_size: usize,
        full_block_size: usize,
        reserved1: *mut c_void,
    }
    type Visitor = unsafe extern "C" fn(heap: *const c_void, area: *const Area, block: *mut c_void, block_size: usize, arg: *mut c_void) -> bool;
    unsafe extern "C" {
        fn mi_heap_main() -> *mut c_void;
        fn mi_heap_visit_blocks(heap: *mut c_void, visit_blocks: bool, visitor: Visitor, arg: *mut c_void) -> bool;
    }
    // block size -> (areas, committed bytes, used bytes)
    type Totals = BTreeMap<usize, (usize, usize, usize)>;
    unsafe extern "C" fn visit(_heap: *const c_void, area: *const Area, _block: *mut c_void, _block_size: usize, arg: *mut c_void) -> bool {
        let totals = unsafe { &mut *(arg as *mut Totals) };
        let area = unsafe { &*area };
        let entry = totals.entry(area.block_size).or_insert((0, 0, 0));
        entry.0 += 1;
        entry.1 += area.committed;
        entry.2 += area.used * area.full_block_size;
        true
    }
    for (name, heap) in [("main heap", unsafe { mi_heap_main() }), ("read heap", crate::assembly::read::read_heap::raw_heap())] {
        if heap.is_null() {
            continue;
        }
        let mut totals: Totals = BTreeMap::new();
        unsafe { mi_heap_visit_blocks(heap, false, visit, &mut totals as *mut Totals as *mut c_void) };
        const MB: f64 = 1024.0 * 1024.0;
        let committed: usize = totals.values().map(|x| x.1).sum();
        let used: usize = totals.values().map(|x| x.2).sum();
        let mut classes: Vec<(usize, (usize, usize, usize))> = totals.into_iter().collect();
        classes.sort_by_key(|x| std::cmp::Reverse(x.1 .1));
        let top: Vec<String> = classes.iter().take(10).map(|(size, (_, c, u))| format!("{}:{:.0}/{:.0}", size, *u as f64 / MB, *c as f64 / MB)).collect();
        crate::info!("memory {} at {}: used({:.0}MB) committed({:.0}MB) blockSize:usedMB/committedMB {}", name, stage, used as f64 / MB, committed as f64 / MB, top.join(" "));
    }
}
#[cfg(not(feature = "mimalloc"))]
pub fn log_allocator_size_classes(_stage: &str) {}
