//! FFI-safe types for Gelfgren library.
//!
//! All types use `#[repr(C)]` for C-compatible memory layout.

use gelfgren_core::bernstein::BernsteinPolynomial;
use gelfgren_core::rational::RationalFunction;
use gelfgren_core::pade::PadeApproximant;
use gelfgren_core::hermite::LagrangeHermiteInterpolant;
use gelfgren_core::piecewise::{Mesh, PiecewiseRational};
use gelfgren_core::bvp::BoundaryValueProblem;

/// Opaque type representing a Bernstein polynomial.
///
/// Create with `gelfgren_bernstein_create()`, free with `gelfgren_bernstein_free()`.
#[repr(C)]
pub struct GelfgrenBernstein {
    _private: [u8; 0],
}

/// Opaque type representing a rational function.
///
/// Create with `gelfgren_rational_create()`, free with `gelfgren_rational_free()`.
#[repr(C)]
pub struct GelfgrenRational {
    _private: [u8; 0],
}

/// Opaque type representing a Padé approximant.
///
/// Create with `gelfgren_pade_from_series()`, free with `gelfgren_pade_free()`.
#[repr(C)]
pub struct GelfgrenPade {
    _private: [u8; 0],
}

/// Opaque type representing a Hermite interpolant.
///
/// Create with `gelfgren_hermite_create()`, free with `gelfgren_hermite_free()`.
#[repr(C)]
pub struct GelfgrenHermite {
    _private: [u8; 0],
}

/// Opaque type representing a mesh partition.
///
/// Create with `gelfgren_mesh_uniform()` or `gelfgren_mesh_chebyshev()`,
/// free with `gelfgren_mesh_free()`.
#[repr(C)]
pub struct GelfgrenMesh {
    _private: [u8; 0],
}

/// Opaque type representing a piecewise rational function.
///
/// Create with `gelfgren_piecewise_create()`, free with `gelfgren_piecewise_free()`.
#[repr(C)]
pub struct GelfgrenPiecewise {
    _private: [u8; 0],
}

/// Opaque type representing a boundary value problem.
///
/// Create with `gelfgren_bvp_create()`, free with `gelfgren_bvp_free()`.
#[repr(C)]
pub struct GelfgrenBVP {
    _private: [u8; 0],
}

/// Error codes returned by Gelfgren FFI functions.
///
/// Zero indicates success, negative values indicate errors.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GelfgrenErrorCode {
    /// Operation succeeded
    Success = 0,
    /// Null pointer provided as argument
    NullPointer = -1,
    /// Invalid dimension or size
    InvalidDimension = -2,
    /// Invalid argument value
    InvalidArgument = -3,
    /// Division by zero (pole in rational function)
    DivisionByZero = -4,
    /// Degree mismatch in operation
    DegreeMismatch = -5,
    /// Interval mismatch in operation
    IntervalMismatch = -6,
    /// Singular matrix (cannot invert)
    SingularMatrix = -7,
    /// Feature not yet implemented
    NotImplemented = -8,
    /// Unknown or unspecified error
    Unknown = -99,
}

/// Boundary condition type for BVPs.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GelfgrenBoundaryType {
    /// Dirichlet: u(x) = value
    Dirichlet = 0,
    /// Neumann: u'(x) = value
    Neumann = 1,
    /// Robin: c₀·u(x) + c₁·u'(x) = value
    Robin = 2,
}

/// Boundary condition specification for BVPs.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct GelfgrenBoundaryCondition {
    /// Type of boundary condition
    pub bc_type: GelfgrenBoundaryType,
    /// Location (x-coordinate)
    pub location: f64,
    /// Value on right-hand side
    pub value: f64,
    /// Coefficients [c₀, c₁] for Robin BC (unused for Dirichlet/Neumann)
    pub coefficients: [f64; 2],
}

// Internal wrapper types to hold actual Rust objects
pub(crate) struct BernsteinBox {
    pub poly: BernsteinPolynomial<f64>,
}

pub(crate) struct RationalBox {
    pub rational: RationalFunction<f64>,
}

pub(crate) struct PadeBox {
    pub pade: PadeApproximant<f64>,
}

pub(crate) struct HermiteBox {
    pub hermite: LagrangeHermiteInterpolant<f64>,
}

pub(crate) struct MeshBox {
    pub mesh: Mesh<f64>,
}

pub(crate) struct PiecewiseBox {
    pub piecewise: PiecewiseRational<f64>,
}

pub(crate) struct BVPBox {
    pub bvp: BoundaryValueProblem<f64>,
}

// Conversion functions between opaque pointers and boxed types
impl GelfgrenBernstein {
    pub(crate) fn from_box(b: Box<BernsteinBox>) -> *mut Self {
        Box::into_raw(b) as *mut Self
    }

    pub(crate) unsafe fn as_box<'a>(ptr: *mut Self) -> &'a mut BernsteinBox {
        &mut *(ptr as *mut BernsteinBox)
    }

    pub(crate) unsafe fn into_box(ptr: *mut Self) -> Box<BernsteinBox> {
        Box::from_raw(ptr as *mut BernsteinBox)
    }
}

impl GelfgrenRational {
    pub(crate) fn from_box(b: Box<RationalBox>) -> *mut Self {
        Box::into_raw(b) as *mut Self
    }

    pub(crate) unsafe fn as_box<'a>(ptr: *mut Self) -> &'a mut RationalBox {
        &mut *(ptr as *mut RationalBox)
    }

    pub(crate) unsafe fn into_box(ptr: *mut Self) -> Box<RationalBox> {
        Box::from_raw(ptr as *mut RationalBox)
    }
}

impl GelfgrenPade {
    pub(crate) fn from_box(b: Box<PadeBox>) -> *mut Self {
        Box::into_raw(b) as *mut Self
    }

    pub(crate) unsafe fn as_box<'a>(ptr: *mut Self) -> &'a mut PadeBox {
        &mut *(ptr as *mut PadeBox)
    }

    pub(crate) unsafe fn into_box(ptr: *mut Self) -> Box<PadeBox> {
        Box::from_raw(ptr as *mut PadeBox)
    }
}

impl GelfgrenHermite {
    pub(crate) fn from_box(b: Box<HermiteBox>) -> *mut Self {
        Box::into_raw(b) as *mut Self
    }

    pub(crate) unsafe fn as_box<'a>(ptr: *mut Self) -> &'a mut HermiteBox {
        &mut *(ptr as *mut HermiteBox)
    }

    pub(crate) unsafe fn into_box(ptr: *mut Self) -> Box<HermiteBox> {
        Box::from_raw(ptr as *mut HermiteBox)
    }
}

impl GelfgrenMesh {
    pub(crate) fn from_box(b: Box<MeshBox>) -> *mut Self {
        Box::into_raw(b) as *mut Self
    }

    pub(crate) unsafe fn as_box<'a>(ptr: *mut Self) -> &'a mut MeshBox {
        &mut *(ptr as *mut MeshBox)
    }

    pub(crate) unsafe fn into_box(ptr: *mut Self) -> Box<MeshBox> {
        Box::from_raw(ptr as *mut MeshBox)
    }
}

impl GelfgrenPiecewise {
    pub(crate) fn from_box(b: Box<PiecewiseBox>) -> *mut Self {
        Box::into_raw(b) as *mut Self
    }

    pub(crate) unsafe fn as_box<'a>(ptr: *mut Self) -> &'a mut PiecewiseBox {
        &mut *(ptr as *mut PiecewiseBox)
    }

    pub(crate) unsafe fn into_box(ptr: *mut Self) -> Box<PiecewiseBox> {
        Box::from_raw(ptr as *mut PiecewiseBox)
    }
}

impl GelfgrenBVP {
    pub(crate) fn from_box(b: Box<BVPBox>) -> *mut Self {
        Box::into_raw(b) as *mut Self
    }

    pub(crate) unsafe fn as_box<'a>(ptr: *mut Self) -> &'a mut BVPBox {
        &mut *(ptr as *mut BVPBox)
    }

    pub(crate) unsafe fn into_box(ptr: *mut Self) -> Box<BVPBox> {
        Box::from_raw(ptr as *mut BVPBox)
    }
}
