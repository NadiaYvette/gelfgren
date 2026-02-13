//! Memory management utilities for FFI.
//!
//! Provides allocation and deallocation functions for all Gelfgren types.

use crate::types::*;

/// Frees a Bernstein polynomial.
///
/// # Safety
///
/// The pointer must have been created by a Gelfgren function and not already freed.
/// After calling this function, the pointer is invalid.
///
/// # Example (C)
///
/// ```c
/// GelfgrenBernstein* poly = gelfgren_bernstein_create(...);
/// // ... use poly ...
/// gelfgren_bernstein_free(poly);
/// ```
#[no_mangle]
pub unsafe extern "C" fn gelfgren_bernstein_free(ptr: *mut GelfgrenBernstein) {
    if !ptr.is_null() {
        let _ = GelfgrenBernstein::into_box(ptr);
        // Box is dropped here, freeing memory
    }
}

/// Frees a rational function.
///
/// # Safety
///
/// The pointer must have been created by a Gelfgren function and not already freed.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_rational_free(ptr: *mut GelfgrenRational) {
    if !ptr.is_null() {
        let _ = GelfgrenRational::into_box(ptr);
    }
}

/// Frees a Padé approximant.
///
/// # Safety
///
/// The pointer must have been created by a Gelfgren function and not already freed.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_pade_free(ptr: *mut GelfgrenPade) {
    if !ptr.is_null() {
        let _ = GelfgrenPade::into_box(ptr);
    }
}

/// Frees a Hermite interpolant.
///
/// # Safety
///
/// The pointer must have been created by a Gelfgren function and not already freed.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_hermite_free(ptr: *mut GelfgrenHermite) {
    if !ptr.is_null() {
        let _ = GelfgrenHermite::into_box(ptr);
    }
}

/// Frees a mesh.
///
/// # Safety
///
/// The pointer must have been created by a Gelfgren function and not already freed.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_mesh_free(ptr: *mut GelfgrenMesh) {
    if !ptr.is_null() {
        let _ = GelfgrenMesh::into_box(ptr);
    }
}

/// Frees a piecewise rational function.
///
/// # Safety
///
/// The pointer must have been created by a Gelfgren function and not already freed.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_piecewise_free(ptr: *mut GelfgrenPiecewise) {
    if !ptr.is_null() {
        let _ = GelfgrenPiecewise::into_box(ptr);
    }
}

/// Frees a boundary value problem.
///
/// # Safety
///
/// The pointer must have been created by a Gelfgren function and not already freed.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_bvp_free(ptr: *mut GelfgrenBVP) {
    if !ptr.is_null() {
        let _ = GelfgrenBVP::into_box(ptr);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gelfgren_core::bernstein::BernsteinPolynomial;

    #[test]
    fn test_bernstein_alloc_free() {
        unsafe {
            // Create a Bernstein polynomial
            let poly = BernsteinPolynomial::new(vec![1.0, 2.0, 3.0], 0.0, 1.0).unwrap();
            let boxed = Box::new(BernsteinBox { poly });
            let ptr = GelfgrenBernstein::from_box(boxed);

            // Verify it's not null
            assert!(!ptr.is_null());

            // Free it
            gelfgren_bernstein_free(ptr);

            // ptr is now invalid (can't test this without UB, but at least no crash)
        }
    }

    #[test]
    fn test_null_free() {
        unsafe {
            // Should not crash
            gelfgren_bernstein_free(std::ptr::null_mut());
            gelfgren_rational_free(std::ptr::null_mut());
            gelfgren_pade_free(std::ptr::null_mut());
            gelfgren_hermite_free(std::ptr::null_mut());
            gelfgren_mesh_free(std::ptr::null_mut());
            gelfgren_piecewise_free(std::ptr::null_mut());
            gelfgren_bvp_free(std::ptr::null_mut());
        }
    }
}
