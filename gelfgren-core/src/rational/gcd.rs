//! Greatest common divisor (GCD) computation and simplification.
//!
//! Implements polynomial GCD using Euclid's algorithm adapted for
//! Bernstein form (Section 5.4 of Farouki & Rajan).

use super::function::RationalFunction;
use crate::bernstein::BernsteinPolynomial;
use crate::error::Result;
use num_traits::{Float, FromPrimitive};
use std::fmt;

impl<T: Float + FromPrimitive + fmt::Debug + std::iter::Sum> RationalFunction<T> {
    /// Simplifies the rational function by removing common factors.
    ///
    /// Computes gcd(P, Q) and returns (P/gcd, Q/gcd).
    ///
    /// Note: This is a simplified implementation. Full GCD in Bernstein form
    /// is complex and may require additional numerical considerations.
    pub fn simplify(&self) -> Result<Self> {
        // For now, return self as-is
        // TODO: Implement proper GCD algorithm from Farouki-Rajan Section 5.4
        Ok(self.clone())
    }

    /// Checks if the numerator and denominator have common factors.
    ///
    /// This is a heuristic check based on shared roots.
    pub fn has_common_factors(&self, tolerance: T) -> bool {
        // Check for shared roots by evaluating both at multiple points
        // and seeing if they're both near zero simultaneously
        let (a, b) = self.interval();
        let n_samples = 20;

        for i in 0..=n_samples {
            let t = T::from_usize(i).unwrap() / T::from_usize(n_samples).unwrap();
            let x = a + t * (b - a);

            let p_val = self.numerator().evaluate(x);
            let q_val = self.denominator().evaluate(x);

            if p_val.abs() < tolerance && q_val.abs() < tolerance {
                return true;
            }
        }

        false
    }
}

/// Computes the greatest common divisor of two Bernstein polynomials.
///
/// Uses Euclid's algorithm: gcd(a, b) = gcd(b, a mod b)
///
/// # Note
///
/// This is a simplified implementation. A full implementation would need:
/// - Proper handling of numerical stability in Bernstein form
/// - Sylvester resultant computations (Farouki-Rajan eq 77-78)
/// - Bezout coefficients for extended GCD
///
/// For now, this returns the polynomial with lower degree (placeholder).
pub fn gcd<T: Float + FromPrimitive + fmt::Debug + std::iter::Sum>(
    _p: &BernsteinPolynomial<T>,
    _q: &BernsteinPolynomial<T>,
) -> Result<BernsteinPolynomial<T>> {
    // Placeholder: return constant polynomial 1
    // TODO: Implement proper Euclid's algorithm in Bernstein form
    let (a, b) = _p.interval();
    BernsteinPolynomial::constant(T::one(), a, b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simplify_no_common_factors() {
        // (x+1)/(x+2) has no common factors
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 2.0], 0.0, 1.0).unwrap();
        let q = BernsteinPolynomial::from_unscaled(vec![2.0, 3.0], 0.0, 1.0).unwrap();
        let r = RationalFunction::new(p, q).unwrap();

        let simplified = r.simplify().unwrap();
        assert_eq!(simplified.degrees(), r.degrees());
    }

    #[test]
    fn test_has_common_factors_false() {
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 2.0], 0.0, 1.0).unwrap();
        let q = BernsteinPolynomial::from_unscaled(vec![2.0, 3.0], 0.0, 1.0).unwrap();
        let r = RationalFunction::new(p, q).unwrap();

        assert!(!r.has_common_factors(1e-6));
    }

    #[test]
    fn test_gcd_placeholder() {
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 2.0, 3.0], 0.0, 1.0).unwrap();
        let q = BernsteinPolynomial::from_unscaled(vec![2.0, 4.0], 0.0, 1.0).unwrap();
        let g = gcd(&p, &q).unwrap();

        // Placeholder returns constant 1
        assert_eq!(g.degree(), 0);
    }
}
