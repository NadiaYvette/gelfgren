//! FFI functions for piecewise rational functions.

use gelfgren_core::piecewise::Mesh;
use crate::types::{GelfgrenMesh, GelfgrenPiecewise, GelfgrenErrorCode, MeshBox};
use crate::error::{result_to_code, set_last_error};

/// Creates a uniform mesh on interval [a, b] with n subintervals.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_mesh_uniform(
    a: f64,
    b: f64,
    n: usize,
) -> *mut GelfgrenMesh {
    let result = Mesh::uniform(a, b, n);
    let (mesh_opt, code) = result_to_code(result);

    if code != GelfgrenErrorCode::Success {
        return std::ptr::null_mut();
    }

    let mesh = mesh_opt.unwrap();
    let boxed = Box::new(MeshBox { mesh });
    GelfgrenMesh::from_box(boxed)
}

/// Creates a Chebyshev mesh on interval [a, b] with n subintervals.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_mesh_chebyshev(
    a: f64,
    b: f64,
    n: usize,
) -> *mut GelfgrenMesh {
    let result = Mesh::chebyshev(a, b, n);
    let (mesh_opt, code) = result_to_code(result);

    if code != GelfgrenErrorCode::Success {
        return std::ptr::null_mut();
    }

    let mesh = mesh_opt.unwrap();
    let boxed = Box::new(MeshBox { mesh });
    GelfgrenMesh::from_box(boxed)
}

/// Evaluates a piecewise rational function at a point.
#[no_mangle]
pub unsafe extern "C" fn gelfgren_piecewise_evaluate(
    piecewise: *const GelfgrenPiecewise,
    x: f64,
    result: *mut f64,
) -> GelfgrenErrorCode {
    if piecewise.is_null() {
        set_last_error("Null pointer: piecewise".to_string());
        return GelfgrenErrorCode::NullPointer;
    }
    if result.is_null() {
        set_last_error("Null pointer: result".to_string());
        return GelfgrenErrorCode::NullPointer;
    }

    let piecewise_box = GelfgrenPiecewise::as_box(piecewise as *mut _);
    let eval_result = piecewise_box.piecewise.evaluate(x);
    let (value_opt, code) = result_to_code(eval_result);

    if code != GelfgrenErrorCode::Success {
        return code;
    }

    *result = value_opt.unwrap();
    GelfgrenErrorCode::Success
}
