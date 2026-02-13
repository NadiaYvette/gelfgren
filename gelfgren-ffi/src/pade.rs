//! FFI functions for Padé approximants.

use std::slice;
use gelfgren_core::pade::PadeApproximant;
use crate::types::{GelfgrenPade, GelfgrenErrorCode, PadeBox};
use crate::error::{result_to_code, set_last_error};

/// Creates a Padé approximant [n/m] from power series coefficients.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_pade_from_series(
    coeffs: *const f64,
    n: usize,
    m: usize,
    center: f64,
    a: f64,
    b: f64,
) -> *mut GelfgrenPade {
    if coeffs.is_null() {
        set_last_error("Null pointer: coeffs".to_string());
        return std::ptr::null_mut();
    }

    let coeffs_slice = slice::from_raw_parts(coeffs, n + m + 1);
    let coeffs_vec = coeffs_slice.to_vec();

    let result = PadeApproximant::from_series(&coeffs_vec, n, m, center, (a, b));
    let (pade_opt, code) = result_to_code(result);

    if code != GelfgrenErrorCode::Success {
        return std::ptr::null_mut();
    }

    let pade = pade_opt.unwrap();
    let boxed = Box::new(PadeBox { pade });
    GelfgrenPade::from_box(boxed)
}

/// Evaluates a Padé approximant at a point.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_pade_evaluate(
    pade: *const GelfgrenPade,
    x: f64,
    result: *mut f64,
) -> GelfgrenErrorCode {
    if pade.is_null() {
        set_last_error("Null pointer: pade".to_string());
        return GelfgrenErrorCode::NullPointer;
    }
    if result.is_null() {
        set_last_error("Null pointer: result".to_string());
        return GelfgrenErrorCode::NullPointer;
    }

    let pade_box = GelfgrenPade::as_box(pade as *mut _);
    let eval_result = pade_box.pade.evaluate(x);
    let (value_opt, code) = result_to_code(eval_result);

    if code != GelfgrenErrorCode::Success {
        return code;
    }

    *result = value_opt.unwrap();
    GelfgrenErrorCode::Success
}
