//! Gelfgren Core - High-performance numerical computing library
//!
//! This library provides numerical computing primitives including:
//! - Linear algebra operations (matrices, vectors, decompositions)
//! - Statistical functions (descriptive statistics, distributions)
//! - Numerical methods (integration, differentiation, ODE solvers)
//! - Special mathematical functions (gamma, bessel, etc.)

pub mod error;
pub mod linalg;
pub mod stats;
pub mod numerical;
pub mod special;

pub use error::{GelfgrenError, Result};
