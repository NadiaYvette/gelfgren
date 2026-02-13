//! Lagrange-Hermite interpolation.
//!
//! Implementation of J.F. Traub's 1964 formulas for interpolating polynomials
//! that match function values and derivatives at multiple points.
//!
//! # Mathematical Background
//!
//! The Lagrange-Hermite problem: Given p(n+1) values y_i^(m) for i=0,...,n
//! and m=0,...,p-1, find polynomial P_{n,p}(t) of degree p(n+1)-1 such that:
//!
//! P_{n,p}^(m)(x_i) = y_i^(m)  for all i,m
//!
//! # Key Components
//!
//! - **Bell Polynomials**: B_n(ω; g_1,...,g_n) for derivative transformations
//! - **Lagrange Basis**: L_i^p(t) standard Lagrange polynomials
//! - **Traub's Formula**: Equation 3.6 combining these elements
//!
//! # Applications
//!
//! - Hermite interpolation (function + first derivative)
//! - Higher-order contact conditions
//! - Padé approximant construction from derivative data
//! - Gelfgren's symmetric approximants with endpoint derivatives
//!
//! # References
//!
//! - Traub (1964): "On Lagrange-Hermite Interpolation"
//! - Bell (1927): Original definition of exponential polynomials
//! - Comtet (1974): "Advanced Combinatorics" for Bell polynomial properties

mod bell;
mod interpolation;

pub use bell::BellPolynomial;
pub use interpolation::{HermiteData, LagrangeHermiteInterpolant};
