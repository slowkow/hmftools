//! Shared mutable references in the manner of Java object references.
//!
//! The Java code is a graph of objects that refer to and change each other (assemblies, phase groups, links,
//! reads and their mates), compared by identity. `P<T>` is a copyable pointer to a value that lives for the
//! whole run; `A<T>` is a counted pointer used for reads so their memory is returned once nothing refers to
//! them. Both hand out `&mut T` from a shared handle through `m()`. As in Java, it is up to the code to keep
//! objects that are being changed on one thread at a time; borrows from `m()` must be short and must not be
//! held across calls that can reach the same object.

use std::cell::UnsafeCell;
use std::hash::{Hash, Hasher};
use std::ops::Deref;
use std::ptr::NonNull;
use std::sync::Arc;

pub struct P<T>(NonNull<T>);

impl<T> P<T> {
    pub fn new(value: T) -> Self {
        P(NonNull::from(Box::leak(Box::new(value))))
    }
    #[allow(clippy::mut_from_ref)]
    #[inline]
    pub fn m(&self) -> &mut T {
        unsafe { &mut *self.0.as_ptr() }
    }
    #[inline]
    pub fn addr(&self) -> usize {
        self.0.as_ptr() as usize
    }
}
impl<T> Clone for P<T> {
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for P<T> {}
impl<T> Deref for P<T> {
    type Target = T;
    #[inline]
    fn deref(&self) -> &T {
        unsafe { &*self.0.as_ptr() }
    }
}
impl<T> PartialEq for P<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}
impl<T> Eq for P<T> {}
impl<T> Hash for P<T> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.addr().hash(state)
    }
}
unsafe impl<T> Send for P<T> {}
unsafe impl<T> Sync for P<T> {}

struct Cell<T>(UnsafeCell<T>);
unsafe impl<T> Send for Cell<T> {}
unsafe impl<T> Sync for Cell<T> {}

pub struct A<T>(Arc<Cell<T>>);

impl<T> A<T> {
    pub fn new(value: T) -> Self {
        A(Arc::new(Cell(UnsafeCell::new(value))))
    }
    #[allow(clippy::mut_from_ref)]
    #[inline]
    pub fn m(&self) -> &mut T {
        unsafe { &mut *self.0.0.get() }
    }
    #[inline]
    pub fn addr(&self) -> usize {
        Arc::as_ptr(&self.0) as usize
    }
    pub fn strong_count(&self) -> usize {
        Arc::strong_count(&self.0)
    }
}
impl<T> Clone for A<T> {
    #[inline]
    fn clone(&self) -> Self {
        A(self.0.clone())
    }
}
impl<T> Deref for A<T> {
    type Target = T;
    #[inline]
    fn deref(&self) -> &T {
        unsafe { &*self.0.0.get() }
    }
}
impl<T> PartialEq for A<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}
impl<T> Eq for A<T> {}
impl<T> Hash for A<T> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.addr().hash(state)
    }
}

/// identity of two optional references (Java `a == b` with nulls)
pub fn same<T: PartialEq>(a: &Option<T>, b: &Option<T>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(x), Some(y)) => x == y,
        _ => false,
    }
}

impl<T> crate::common::java_id_hash::IdKey for P<T> {
    fn id_addr(&self) -> usize {
        self.addr()
    }
}
impl<T> crate::common::java_id_hash::IdKey for A<T> {
    fn id_addr(&self) -> usize {
        self.addr()
    }
}
