//! FFI functions for Hermite interpolation.

// Placeholder implementation
use crate::types::{GelfgrenHermite, GelfgrenErrorCode};
use crate::{check_null, ffi_catch};

/// Evaluates a Hermite interpolant at a point.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_hermite_evaluate(
    hermite: *const GelfgrenHermite,
    _x: f64,
    result: *mut f64,
) -> GelfgrenErrorCode {
    {
        check_null!(hermite, "hermite");
        check_null!(result, "result");

        // TODO: Implement Hermite interpolation FFI
        *result = 0.0;
        GelfgrenErrorCode::NotImplemented
    }
}
