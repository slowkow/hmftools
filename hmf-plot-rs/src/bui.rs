//! The BUI software raster and text engine, the C API of `gpu/include/GPU_raster.h` in BUI's kit.
//!
//! `blenderui-sys` builds and links the kit (static libraries bui_kit, bui_gpu, bui_core) and declares part
//! of this API; the rest is declared here against the same libraries. Only this module and `raster.rs`
//! touch BUI: everything else in the crate draws into a display list.

#![allow(non_camel_case_types, dead_code)]

use std::ffi::{c_char, c_float, c_int, c_void};

pub use blenderui_sys::BuiRaster;

pub const BUI_FONT_BOLD: c_int = 1;
pub const BUI_FONT_ITALIC: c_int = 2;

pub const BUI_MESH_INDEX32: c_int = 1;
pub const BUI_MESH_NEAREST: c_int = 2;

/// One textured mesh vertex (the layout of Dear ImGui's ImDrawVert): `col` is 0xAABBGGRR, straight alpha.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct BuiMeshVertex {
    pub x: c_float,
    pub y: c_float,
    pub u: c_float,
    pub v: c_float,
    pub col: u32,
}

#[repr(C)]
pub struct BuiTextRun {
    _private: [u8; 0],
}

unsafe extern "C" {
    pub fn BUI_raster_create(width: c_int, height: c_int) -> *mut BuiRaster;
    pub fn BUI_raster_free(r: *mut BuiRaster);
    pub fn BUI_raster_begin_frame(r: *mut BuiRaster, width: c_int, height: c_int);
    pub fn BUI_raster_end_frame(r: *mut BuiRaster);
    pub fn BUI_raster_width(r: *const BuiRaster) -> c_int;
    pub fn BUI_raster_height(r: *const BuiRaster) -> c_int;
    pub fn BUI_raster_pixels(r: *const BuiRaster) -> *const u32;
    pub fn BUI_raster_set_clip(r: *mut BuiRaster, x: c_int, y: c_int, w: c_int, h: c_int);
    pub fn BUI_raster_fill_rect(r: *mut BuiRaster, x: c_float, y: c_float, w: c_float, h: c_float, rgba: *const c_float);
    pub fn BUI_raster_circle(r: *mut BuiRaster, cx: c_float, cy: c_float, radius: c_float, line_width: c_float, rgba: *const c_float);
    pub fn BUI_raster_line(r: *mut BuiRaster, x1: c_float, y1: c_float, x2: c_float, y2: c_float, width: c_float, rgba: *const c_float);
    pub fn BUI_raster_polygon(r: *mut BuiRaster, points: *const [c_float; 2], n: c_int, rgba: *const c_float);
    pub fn BUI_raster_mesh(
        r: *mut BuiRaster,
        vtx: *const BuiMeshVertex,
        n_vtx: c_int,
        idx: *const c_void,
        n_idx: c_int,
        flags: c_int,
        tex: c_int,
        clip: *const c_float,
        ox: c_float,
        oy: c_float,
    );
    pub fn BUI_raster_set_font(r: *mut BuiRaster, family: *const c_char, flags: c_int);
    pub fn BUI_raster_text(r: *mut BuiRaster, x: c_float, y: c_float, utf8: *const c_char, len: c_int, size: c_float, rgba: *const c_float) -> c_int;

    pub fn BUI_text_measure(family: *const c_char, flags: c_int, size: c_float, utf8: *const c_char, len: c_int, r_width: *mut c_float) -> c_int;
    pub fn BUI_text_add_font(family: *const c_char, flags: c_int, data: *const c_void, len: c_int) -> c_int;
    pub fn BUI_text_vmetrics(family: *const c_char, flags: c_int, r_ascent: *mut c_float, r_descent: *mut c_float, r_gap: *mut c_float) -> c_int;
    pub fn BUI_text_prepare(family: *const c_char, flags: c_int, size: c_float, utf8: *const c_char, len: c_int, r_width: *mut c_float) -> *const BuiTextRun;
    pub fn BUI_text_run_mask(run: *const BuiTextRun, r_w: *mut c_int, r_h: *mut c_int, r_ox: *mut c_int, r_ascent: *mut c_int) -> *const u8;
    pub fn BUI_text_generation() -> u32;

    pub fn BUI_texture_create(w: c_int, h: c_int, rgba: *const c_void, pitch: c_int) -> c_int;
    pub fn BUI_texture_free(tex: c_int);
}
