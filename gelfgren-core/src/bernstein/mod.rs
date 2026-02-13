//! Bernstein polynomial representation and operations.
//!
//! This module implements algorithms from Farouki & Rajan (1987)
//! "Algorithms for Polynomials in Bernstein Form" for numerically
//! stable polynomial operations.
//!
//! # Key Concept
//!
//! A polynomial of degree n on [a,b] is represented as:
//!
//! P(x) = Σ_{k=0}^n C_k^n b_k^n(t)
//!
//! where t = (x-a)/(b-a) maps [a,b] to [0,1] and
//! b_k^n(t) = (n choose k) t^k (1-t)^{n-k} is the Bernstein basis.
//!
//! We use **scaled coefficients** Ĉ_k^n = (n choose k) C_k^n
//! for improved numerical properties and simpler formulas.

mod polynomial;
mod arithmetic;
mod calculus;

pub use polynomial::{binomial, BernsteinPolynomial};
