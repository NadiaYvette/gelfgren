//! Gelfgren FFI - C Foreign Function Interface for Gelfgren
//!
//! This library provides C-compatible bindings to the Gelfgren numerical
//! computing library. It exposes a C API for:
//!
//! - Bernstein polynomials
//! - Rational functions
//! - Padé approximants
//! - Lagrange-Hermite interpolation
//! - Piecewise rational approximation
//! - Boundary value problem solving
//!
//! # Safety
//!
//! All public FFI functions are `unsafe` from Rust's perspective but designed
//! to be safe when called correctly from C. The API follows these conventions:
//!
//! - All functions return error codes (0 = success, negative = error)
//! - Detailed error messages available via `gelfgren_last_error_message()`
//! - Memory management through explicit create/free pairs
//! - Null pointer checking on all inputs
//! - Thread-local error storage
//!
//! # Error Handling
//!
//! ```c
//! GelfgrenBernstein* poly = gelfgren_bernstein_create(coeffs, 3, 0.0, 1.0);
//! if (poly == NULL) {
//!     const char* error = gelfgren_last_error_message();
//!     fprintf(stderr, "Error: %s\n", error);
//!     return -1;
//! }
//! // ... use poly ...
//! gelfgren_bernstein_free(poly);
//! ```

mod types;
mod error;
mod memory;
mod bernstein;
mod rational;
mod pade;
mod hermite;
mod piecewise;
mod bvp;

pub use types::*;
pub use error::*;
pub use memory::*;
pub use bernstein::*;
pub use rational::*;
pub use pade::*;
pub use hermite::*;
pub use piecewise::*;
pub use bvp::*;
