//! A separate allocator heap for reads (not part of the Java code).
//!
//! Nearly all memory the step allocates is for reads that are freed again soon after. When they come from
//! the same heap as everything else, the few long-lived allocations of the same sizes made in between keep
//! the allocator's pages from emptying, and the freed memory is not returned to the system. The read
//! objects and the buffers of their records are therefore allocated here, from a heap that holds nothing
//! else, whose pages empty when the reads are freed.
//!
//! This uses the heap interface of the mimalloc allocator: a block allocated from any of its heaps is freed
//! by the ordinary free function, which is what the global allocator of this program calls, so the buffers
//! can be owned by ordinary `Vec` and `String` values. With another global allocator (the `mimalloc`
//! feature turned off) everything comes from the global allocator.

use std::alloc::Layout;

#[cfg(feature = "mimalloc")]
mod mi {
    use std::ffi::c_void;
    use std::sync::OnceLock;

    unsafe extern "C" {
        fn mi_version() -> i32;
        fn mi_heap_new() -> *mut c_void;
        fn mi_heap_malloc(heap: *mut c_void, size: usize) -> *mut c_void;
        fn mi_heap_malloc_aligned(heap: *mut c_void, size: usize, alignment: usize) -> *mut c_void;
        fn mi_heap_collect(heap: *mut c_void, force: bool);
    }

    struct Heap(*mut c_void);
    unsafe impl Send for Heap {}
    unsafe impl Sync for Heap {}
    static READ_HEAP: OnceLock<Heap> = OnceLock::new();

    /// the heap if it has been created
    pub fn heap() -> *mut c_void {
        READ_HEAP.get().map_or(std::ptr::null_mut(), |h| h.0)
    }

    pub fn collect() {
        let heap = heap();
        if !heap.is_null() {
            unsafe { mi_heap_collect(heap, true) };
        }
    }

    pub fn alloc(size: usize, alignment: usize) -> *mut u8 {
        // only version 3 heaps may be used from every thread (a version 2 heap belongs to the thread that
        // created it); with an older allocator the reads come from the global allocator
        let heap = READ_HEAP.get_or_init(|| {
            Heap(if unsafe { mi_version() } >= 30000 { unsafe { mi_heap_new() } } else { std::ptr::null_mut() })
        });
        if heap.0.is_null() {
            return std::ptr::null_mut();
        }
        // every block is aligned to 8 bytes, and to 16 from that size on; asking for an alignment makes the
        // allocator reserve extra room
        if alignment <= 8 || (alignment <= 16 && size >= 16) {
            unsafe { mi_heap_malloc(heap.0, size) as *mut u8 }
        } else {
            unsafe { mi_heap_malloc_aligned(heap.0, size, alignment) as *mut u8 }
        }
    }
}

/// Returns the free memory of the read heap to the system (called when most reads have been freed).
pub fn release_free_memory() {
    #[cfg(feature = "mimalloc")]
    mi::collect();
}

/// The allocator heap that holds the reads, or null (for the memory diagnostics).
#[cfg(feature = "mimalloc")]
pub fn raw_heap() -> *mut std::ffi::c_void {
    mi::heap()
}

/// Allocates from the read heap; the block is freed with `dealloc` or by the global allocator.
pub fn alloc(layout: Layout) -> *mut u8 {
    debug_assert!(layout.size() > 0);
    #[cfg(feature = "mimalloc")]
    {
        let p = mi::alloc(layout.size(), layout.align());
        if !p.is_null() {
            return p;
        }
    }
    let p = unsafe { std::alloc::alloc(layout) };
    if p.is_null() {
        std::alloc::handle_alloc_error(layout);
    }
    p
}

/// Frees a block from `alloc`.
///
/// # Safety
/// `p` must come from `alloc` with the same layout and must not be used afterwards.
pub unsafe fn dealloc(p: *mut u8, layout: Layout) {
    unsafe { std::alloc::dealloc(p, layout) }
}

/// A copy of `source` in the read heap.
pub fn vec_from<T: Copy>(source: &[T]) -> Vec<T> {
    if source.is_empty() || std::mem::size_of::<T>() == 0 {
        return Vec::new();
    }
    let layout = Layout::array::<T>(source.len()).expect("layout");
    let p = alloc(layout) as *mut T;
    unsafe {
        std::ptr::copy_nonoverlapping(source.as_ptr(), p, source.len());
        Vec::from_raw_parts(p, source.len(), source.len())
    }
}

/// A copy of `source` in the read heap.
pub fn string_from(source: &str) -> String {
    unsafe { String::from_utf8_unchecked(vec_from(source.as_bytes())) }
}
