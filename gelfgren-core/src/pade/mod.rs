//! Padé approximant construction.
//!
//! A Padé approximant [n/m] is a rational function P_n(x)/Q_m(x) that
//! matches a given power series (or function with derivatives) to high order.
//!
//! # Mathematical Background
//!
//! Given a function f(x) with power series expansion:
//! f(x) = a₀ + a₁x + a₂x² + ... + aₖxᵏ + ...
//!
//! The [n/m] Padé approximant R(x) = Pₙ(x)/Qₘ(x) satisfies:
//! f(x) - R(x) = O(x^(n+m+1))
//!
//! This means R(x) matches f(x) and its first n+m derivatives at x=0.
//!
//! # Construction Methods
//!
//! 1. **From power series**: Given coefficients {a₀, a₁, ..., aₙ₊ₘ}
//! 2. **From derivatives**: Given f(x₀), f'(x₀), ..., f^(n+m)(x₀)
//! 3. **Symmetric**: Two-point approximation at x₀ and x₁ (Gelfgren's method)
//!
//! # References
//!
//! - Gelfgren (1975): Symmetric Padé approximants on intervals
//! - Baker & Graves-Morris: "Padé Approximants" (standard reference)

mod approximant;
mod symmetric;

pub use approximant::PadeApproximant;
pub use symmetric::SymmetricPade;
