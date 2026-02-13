//! Piecewise rational function construction.
//!
//! Implementation of Gelfgren's 1975 algorithm for constructing piecewise
//! rational approximations S_{n,m}(τ)(z) on a mesh partition τ.
//!
//! # Mathematical Background
//!
//! Given a partition τ: a = x₀ < x₁ < ... < x_N = b, construct on each
//! subinterval Δⱼ = [xⱼ₋₁, xⱼ] a rational approximant Rⱼ(x) = P_{n,j}(x)/Q_{m,j}(x).
//!
//! ## Key Properties
//!
//! - Each Rⱼ is a symmetric Padé approximant matching derivatives at endpoints
//! - If gcd(P_{n,j}, Q_{m,j}) = 1, then R has C^(p-1) continuity at nodes
//! - Global approximation error depends on mesh size h = max|Δⱼ|
//! - Convergence requires n ≥ k·m for some k > 0
//!
//! ## Mesh Refinement
//!
//! For better approximation:
//! - Uniform refinement: h → h/2
//! - Adaptive refinement: subdivide intervals with large error
//! - Graded meshes: finer spacing near singularities
//!
//! # Algorithm Overview
//!
//! 1. Create mesh partition τ
//! 2. For each subinterval Δⱼ:
//!    a. Collect function/derivative data at endpoints
//!    b. Construct symmetric Padé approximant
//!    c. Check for common factors (simplify if needed)
//! 3. Assemble piecewise function
//! 4. Verify continuity at mesh points
//!
//! # References
//!
//! - Gelfgren (1975): "Piecewise Rational Interpolation"
//! - de Boor (1978): "A Practical Guide to Splines" (for mesh concepts)

mod mesh;
mod construction;
mod evaluation;

pub use mesh::{Mesh, MeshPoint};
pub use construction::PiecewiseRational;
pub use evaluation::SubintervalData;
