//! Rational function representation and operations.
//!
//! A rational function is a ratio of two polynomials P(x)/Q(x).
//! Both numerator and denominator are represented as Bernstein polynomials
//! for numerical stability.
//!
//! # Key Operations
//!
//! - Evaluation with pole detection
//! - Differentiation using quotient rule
//! - Simplification by removing common factors (GCD)
//! - Arithmetic operations (addition, multiplication)

mod function;
mod arithmetic;
mod gcd;

pub use function::RationalFunction;
