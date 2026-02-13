//! Gelfgren Core - Piecewise Rational Interpolation Library
//!
//! Implementation of Jan Gelfgren's 1975 algorithm for piecewise rational interpolation,
//! using Bernstein polynomial representation (Farouki & Rajan 1987) and
//! Lagrange-Hermite interpolation (Traub 1964).
//!
//! # Core Components
//!
//! - **Bernstein Polynomials**: Numerically stable polynomial representation
//! - **Rational Functions**: P_n(x)/Q_m(x) with automatic simplification
//! - **Padé Approximants**: From power series or derivative data
//! - **Lagrange-Hermite Interpolation**: Using Traub's formulas with Bell polynomials
//! - **Piecewise Rational Approximation**: Gelfgren's method for mesh-based construction
//! - **Boundary Value Problems**: Linear constraint systems for ODE solutions

pub mod error;
pub mod bernstein;
pub mod rational;
pub mod pade;
pub mod hermite;
pub mod piecewise;
pub mod bvp;

pub use error::{GelfgrenError, Result};
