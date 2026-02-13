//! FFI functions for rational functions.

use gelfgren_core::rational::RationalFunction;
use crate::types::{GelfgrenRational, GelfgrenBernstein, GelfgrenErrorCode, RationalBox, BernsteinBox};
use crate::error::{result_to_code, set_last_error};

/// Creates a new rational function from numerator and denominator polynomials.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_rational_create(
    numerator: *const GelfgrenBernstein,
    denominator: *const GelfgrenBernstein,
) -> *mut GelfgrenRational {
    if numerator.is_null() {
        set_last_error("Null pointer: numerator".to_string());
        return std::ptr::null_mut();
    }
    if denominator.is_null() {
        set_last_error("Null pointer: denominator".to_string());
        return std::ptr::null_mut();
    }

    let num_box = GelfgrenBernstein::as_box(numerator as *mut _);
    let den_box = GelfgrenBernstein::as_box(denominator as *mut _);

    let result = RationalFunction::new(
        num_box.poly.clone(),
        den_box.poly.clone(),
    );
    let (rational_opt, code) = result_to_code(result);

    if code != GelfgrenErrorCode::Success {
        return std::ptr::null_mut();
    }

    let rational = rational_opt.unwrap();
    let boxed = Box::new(RationalBox { rational });
    GelfgrenRational::from_box(boxed)
}

/// Evaluates a rational function at a point.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_rational_evaluate(
    rational: *const GelfgrenRational,
    x: f64,
    result: *mut f64,
) -> GelfgrenErrorCode {
    if rational.is_null() {
        set_last_error("Null pointer: rational".to_string());
        return GelfgrenErrorCode::NullPointer;
    }
    if result.is_null() {
        set_last_error("Null pointer: result".to_string());
        return GelfgrenErrorCode::NullPointer;
    }

    let rational_box = GelfgrenRational::as_box(rational as *mut _);
    let eval_result = rational_box.rational.evaluate(x);
    let (value_opt, code) = result_to_code(eval_result);

    if code != GelfgrenErrorCode::Success {
        return code;
    }

    *result = value_opt.unwrap();
    GelfgrenErrorCode::Success
}

/// Computes the derivative of a rational function.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_rational_derivative(
    rational: *const GelfgrenRational,
) -> *mut GelfgrenRational {
    if rational.is_null() {
        set_last_error("Null pointer: rational".to_string());
        return std::ptr::null_mut();
    }

    let rational_box = GelfgrenRational::as_box(rational as *mut _);
    let deriv_result = rational_box.rational.derivative();
    let (deriv_opt, code) = result_to_code(deriv_result);

    if code != GelfgrenErrorCode::Success {
        return std::ptr::null_mut();
    }

    let deriv = deriv_opt.unwrap();
    let boxed = Box::new(RationalBox { rational: deriv });
    GelfgrenRational::from_box(boxed)
}
