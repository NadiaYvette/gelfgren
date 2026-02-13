//! FFI functions for Bernstein polynomials.

use std::slice;
use gelfgren_core::bernstein::BernsteinPolynomial;
use crate::types::{GelfgrenBernstein, GelfgrenErrorCode, BernsteinBox};
use crate::error::{result_to_code, set_last_error};

/// Creates a new Bernstein polynomial from unscaled coefficients.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_bernstein_create(
    coeffs: *const f64,
    degree: usize,
    a: f64,
    b: f64,
) -> *mut GelfgrenBernstein {
    if coeffs.is_null() {
        set_last_error("Null pointer: coeffs".to_string());
        return std::ptr::null_mut();
    }

    let n = degree;
    let coeffs_slice = slice::from_raw_parts(coeffs, n + 1);
    let coeffs_vec = coeffs_slice.to_vec();

    let result = BernsteinPolynomial::new(coeffs_vec, a, b);
    let (poly_opt, code) = result_to_code(result);

    if code != GelfgrenErrorCode::Success {
        return std::ptr::null_mut();
    }

    let poly = poly_opt.unwrap();
    let boxed = Box::new(BernsteinBox { poly });
    GelfgrenBernstein::from_box(boxed)
}

/// Evaluates a Bernstein polynomial at a point.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_bernstein_evaluate(
    poly: *const GelfgrenBernstein,
    x: f64,
    result: *mut f64,
) -> GelfgrenErrorCode {
    if poly.is_null() {
        set_last_error("Null pointer: poly".to_string());
        return GelfgrenErrorCode::NullPointer;
    }
    if result.is_null() {
        set_last_error("Null pointer: result".to_string());
        return GelfgrenErrorCode::NullPointer;
    }

    let poly_box = GelfgrenBernstein::as_box(poly as *mut _);
    let value = poly_box.poly.evaluate(x);
    *result = value;

    GelfgrenErrorCode::Success
}

/// Returns the degree of a Bernstein polynomial.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_bernstein_degree(
    poly: *const GelfgrenBernstein,
    result: *mut usize,
) -> GelfgrenErrorCode {
    if poly.is_null() {
        set_last_error("Null pointer: poly".to_string());
        return GelfgrenErrorCode::NullPointer;
    }
    if result.is_null() {
        set_last_error("Null pointer: result".to_string());
        return GelfgrenErrorCode::NullPointer;
    }

    let poly_box = GelfgrenBernstein::as_box(poly as *mut _);
    *result = poly_box.poly.degree();

    GelfgrenErrorCode::Success
}

/// Returns the interval [a, b] for a Bernstein polynomial.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_bernstein_interval(
    poly: *const GelfgrenBernstein,
    a: *mut f64,
    b: *mut f64,
) -> GelfgrenErrorCode {
    if poly.is_null() {
        set_last_error("Null pointer: poly".to_string());
        return GelfgrenErrorCode::NullPointer;
    }
    if a.is_null() {
        set_last_error("Null pointer: a".to_string());
        return GelfgrenErrorCode::NullPointer;
    }
    if b.is_null() {
        set_last_error("Null pointer: b".to_string());
        return GelfgrenErrorCode::NullPointer;
    }

    let poly_box = GelfgrenBernstein::as_box(poly as *mut _);
    let (a_val, b_val) = poly_box.poly.interval();
    *a = a_val;
    *b = b_val;

    GelfgrenErrorCode::Success
}

/// Adds two Bernstein polynomials.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_bernstein_add(
    p1: *const GelfgrenBernstein,
    p2: *const GelfgrenBernstein,
) -> *mut GelfgrenBernstein {
    if p1.is_null() {
        set_last_error("Null pointer: p1".to_string());
        return std::ptr::null_mut();
    }
    if p2.is_null() {
        set_last_error("Null pointer: p2".to_string());
        return std::ptr::null_mut();
    }

    let p1_box = GelfgrenBernstein::as_box(p1 as *mut _);
    let p2_box = GelfgrenBernstein::as_box(p2 as *mut _);

    let add_result = p1_box.poly.clone() + p2_box.poly.clone();
    let (poly_opt, code) = result_to_code(add_result);

    if code != GelfgrenErrorCode::Success {
        return std::ptr::null_mut();
    }

    let poly = poly_opt.unwrap();
    let boxed = Box::new(BernsteinBox { poly });
    GelfgrenBernstein::from_box(boxed)
}

/// Subtracts two Bernstein polynomials.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_bernstein_subtract(
    p1: *const GelfgrenBernstein,
    p2: *const GelfgrenBernstein,
) -> *mut GelfgrenBernstein {
    if p1.is_null() {
        set_last_error("Null pointer: p1".to_string());
        return std::ptr::null_mut();
    }
    if p2.is_null() {
        set_last_error("Null pointer: p2".to_string());
        return std::ptr::null_mut();
    }

    let p1_box = GelfgrenBernstein::as_box(p1 as *mut _);
    let p2_box = GelfgrenBernstein::as_box(p2 as *mut _);

    let sub_result = p1_box.poly.clone() - p2_box.poly.clone();
    let (poly_opt, code) = result_to_code(sub_result);

    if code != GelfgrenErrorCode::Success {
        return std::ptr::null_mut();
    }

    let poly = poly_opt.unwrap();
    let boxed = Box::new(BernsteinBox { poly });
    GelfgrenBernstein::from_box(boxed)
}

/// Multiplies two Bernstein polynomials.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_bernstein_multiply(
    p1: *const GelfgrenBernstein,
    p2: *const GelfgrenBernstein,
) -> *mut GelfgrenBernstein {
    if p1.is_null() {
        set_last_error("Null pointer: p1".to_string());
        return std::ptr::null_mut();
    }
    if p2.is_null() {
        set_last_error("Null pointer: p2".to_string());
        return std::ptr::null_mut();
    }

    let p1_box = GelfgrenBernstein::as_box(p1 as *mut _);
    let p2_box = GelfgrenBernstein::as_box(p2 as *mut _);

    let mul_result = p1_box.poly.clone() * p2_box.poly.clone();
    let (poly_opt, code) = result_to_code(mul_result);

    if code != GelfgrenErrorCode::Success {
        return std::ptr::null_mut();
    }

    let poly = poly_opt.unwrap();
    let boxed = Box::new(BernsteinBox { poly });
    GelfgrenBernstein::from_box(boxed)
}

/// Computes the derivative of a Bernstein polynomial.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_bernstein_derivative(
    poly: *const GelfgrenBernstein,
) -> *mut GelfgrenBernstein {
    if poly.is_null() {
        set_last_error("Null pointer: poly".to_string());
        return std::ptr::null_mut();
    }

    let poly_box = GelfgrenBernstein::as_box(poly as *mut _);
    let deriv_result = poly_box.poly.derivative();
    let (deriv_opt, code) = result_to_code(deriv_result);

    if code != GelfgrenErrorCode::Success {
        return std::ptr::null_mut();
    }

    let deriv = deriv_opt.unwrap();
    let boxed = Box::new(BernsteinBox { poly: deriv });
    GelfgrenBernstein::from_box(boxed)
}

/// Computes the integral (antiderivative) of a Bernstein polynomial.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_bernstein_integral(
    poly: *const GelfgrenBernstein,
) -> *mut GelfgrenBernstein {
    if poly.is_null() {
        set_last_error("Null pointer: poly".to_string());
        return std::ptr::null_mut();
    }

    let poly_box = GelfgrenBernstein::as_box(poly as *mut _);
    let anti_result = poly_box.poly.integral();
    let (anti_opt, code) = result_to_code(anti_result);

    if code != GelfgrenErrorCode::Success {
        return std::ptr::null_mut();
    }

    let anti = anti_opt.unwrap();
    let boxed = Box::new(BernsteinBox { poly: anti });
    GelfgrenBernstein::from_box(boxed)
}

/// Elevates the degree of a Bernstein polynomial by 1.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_bernstein_elevate(
    poly: *const GelfgrenBernstein,
) -> *mut GelfgrenBernstein {
    if poly.is_null() {
        set_last_error("Null pointer: poly".to_string());
        return std::ptr::null_mut();
    }

    let poly_box = GelfgrenBernstein::as_box(poly as *mut _);
    let elevated = poly_box.poly.degree_elevate();
    let boxed = Box::new(BernsteinBox { poly: elevated });
    GelfgrenBernstein::from_box(boxed)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_and_free() {
        unsafe {
            let coeffs = vec![1.0, 2.0, 3.0];
            let poly = gelfgren_bernstein_create(coeffs.as_ptr(), 2, 0.0, 1.0);
            assert!(!poly.is_null());
            crate::memory::gelfgren_bernstein_free(poly);
        }
    }

    #[test]
    fn test_evaluate() {
        unsafe {
            let coeffs = vec![1.0, 2.0, 3.0];
            let poly = gelfgren_bernstein_create(coeffs.as_ptr(), 2, 0.0, 1.0);
            assert!(!poly.is_null());

            let mut result = 0.0;
            let code = gelfgren_bernstein_evaluate(poly, 0.5, &mut result);
            assert_eq!(code, GelfgrenErrorCode::Success);
            assert!(result > 0.0);

            crate::memory::gelfgren_bernstein_free(poly);
        }
    }
}
