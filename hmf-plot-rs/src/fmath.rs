//! The transcendental functions the drawing uses, from the `libm` crate (musl's algorithms in Rust) rather than
//! the platform's C library, so that every platform computes the same bits and writes the same PNG bytes
//! (macOS's and glibc's sin, cos and pow differ in the last bit for some arguments).

pub trait FMath {
    fn fsin(self) -> f64;
    fn fcos(self) -> f64;
    fn fpow(self, e: f64) -> f64;
    fn fpowi(self, e: i32) -> f64;
    fn flog10(self) -> f64;
    fn fcbrt(self) -> f64;
}

impl FMath for f64 {
    #[inline]
    fn fsin(self) -> f64 {
        libm::sin(self)
    }
    #[inline]
    fn fcos(self) -> f64 {
        libm::cos(self)
    }
    #[inline]
    fn fpow(self, e: f64) -> f64 {
        libm::pow(self, e)
    }
    #[inline]
    fn fpowi(self, e: i32) -> f64 {
        libm::pow(self, e as f64)
    }
    #[inline]
    fn flog10(self) -> f64 {
        libm::log10(self)
    }
    #[inline]
    fn fcbrt(self) -> f64 {
        libm::cbrt(self)
    }
}
