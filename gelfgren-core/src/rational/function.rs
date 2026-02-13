//! Core rational function type.

use crate::bernstein::BernsteinPolynomial;
use crate::error::{GelfgrenError, Result};
use num_traits::{Float, FromPrimitive};
use std::fmt;

/// A rational function P(x)/Q(x) where P and Q are Bernstein polynomials.
///
/// # Invariants
///
/// - Numerator and denominator are on the same interval [a, b]
/// - Denominator is not the zero polynomial
/// - Common factors are removed (maintained by simplify())
#[derive(Debug, Clone, PartialEq)]
pub struct RationalFunction<T: Float> {
    /// Numerator polynomial P(x)
    numerator: BernsteinPolynomial<T>,
    /// Denominator polynomial Q(x)
    denominator: BernsteinPolynomial<T>,
}

impl<T: Float + FromPrimitive + fmt::Debug + std::iter::Sum> RationalFunction<T> {
    /// Creates a new rational function P(x)/Q(x).
    ///
    /// # Arguments
    ///
    /// * `numerator` - Numerator polynomial P(x)
    /// * `denominator` - Denominator polynomial Q(x)
    ///
    /// # Errors
    ///
    /// Returns error if:
    /// - Polynomials are on different intervals
    /// - Denominator is the zero polynomial
    pub fn new(
        numerator: BernsteinPolynomial<T>,
        denominator: BernsteinPolynomial<T>,
    ) -> Result<Self> {
        if numerator.interval() != denominator.interval() {
            return Err(GelfgrenError::InvalidArgument(
                "Numerator and denominator must be on the same interval".to_string(),
            ));
        }

        let tol = T::from_f64(1e-10).unwrap();
        if denominator.is_zero(tol) {
            return Err(GelfgrenError::DivisionByZero);
        }

        Ok(Self {
            numerator,
            denominator,
        })
    }

    /// Returns a reference to the numerator polynomial.
    #[inline]
    pub fn numerator(&self) -> &BernsteinPolynomial<T> {
        &self.numerator
    }

    /// Returns a reference to the denominator polynomial.
    #[inline]
    pub fn denominator(&self) -> &BernsteinPolynomial<T> {
        &self.denominator
    }

    /// Returns the interval [a, b] on which the rational function is defined.
    #[inline]
    pub fn interval(&self) -> (T, T) {
        self.numerator.interval()
    }

    /// Returns the degrees (n, m) where numerator has degree n and denominator has degree m.
    #[inline]
    pub fn degrees(&self) -> (usize, usize) {
        (self.numerator.degree(), self.denominator.degree())
    }

    /// Evaluates the rational function at x.
    ///
    /// # Arguments
    ///
    /// * `x` - Point at which to evaluate (must be in [a, b])
    ///
    /// # Returns
    ///
    /// P(x)/Q(x) if Q(x) ≠ 0, otherwise returns an error
    ///
    /// # Errors
    ///
    /// Returns `DivisionByZero` if Q(x) is too close to zero (pole or near-pole)
    pub fn evaluate(&self, x: T) -> Result<T> {
        let p_x = self.numerator.evaluate(x);
        let q_x = self.denominator.evaluate(x);

        let tol = T::from_f64(1e-10).unwrap();
        if q_x.abs() < tol {
            return Err(GelfgrenError::DivisionByZero);
        }

        Ok(p_x / q_x)
    }

    /// Evaluates the rational function at x, returning None if near a pole.
    ///
    /// This is a safe variant that returns None instead of an error for poles.
    pub fn try_evaluate(&self, x: T) -> Option<T> {
        self.evaluate(x).ok()
    }

    /// Checks if x is near a pole of the rational function.
    ///
    /// Returns true if |Q(x)| < tolerance.
    pub fn is_near_pole(&self, x: T, tolerance: T) -> bool {
        self.denominator.evaluate(x).abs() < tolerance
    }

    /// Computes the derivative dR/dx using the quotient rule.
    ///
    /// If R(x) = P(x)/Q(x), then:
    /// dR/dx = (P'(x)Q(x) - P(x)Q'(x)) / Q(x)²
    pub fn derivative(&self) -> Result<Self> {
        let p_prime = self.numerator.derivative()?;
        let q_prime = self.denominator.derivative()?;

        // P'(x)Q(x)
        let mul_result1 = p_prime * self.denominator.clone();
        let term1 = mul_result1?;

        // P(x)Q'(x)
        let mul_result2 = self.numerator.clone() * q_prime;
        let term2 = mul_result2?;

        // P'(x)Q(x) - P(x)Q'(x)
        let sub_result = term1 - term2;
        let new_numerator = sub_result?;

        // Q(x)²
        let mul_result3 = self.denominator.clone() * self.denominator.clone();
        let new_denominator = mul_result3?;

        Self::new(new_numerator, new_denominator)
    }

    /// Computes the k-th derivative.
    pub fn derivative_n(&self, order: usize) -> Result<Self> {
        if order == 0 {
            return Ok(self.clone());
        }

        let mut result = self.derivative()?;
        for _ in 1..order {
            result = result.derivative()?;
        }
        Ok(result)
    }

    /// Creates a polynomial rational function (denominator = 1).
    pub fn from_polynomial(poly: BernsteinPolynomial<T>) -> Result<Self> {
        let (a, b) = poly.interval();
        let one = BernsteinPolynomial::constant(T::one(), a, b)?;
        Self::new(poly, one)
    }

    /// Creates a constant rational function.
    pub fn constant(value: T, a: T, b: T) -> Result<Self> {
        let poly = BernsteinPolynomial::constant(value, a, b)?;
        Self::from_polynomial(poly)
    }

    /// Checks if the rational function is numerically a polynomial (denominator ≈ constant).
    pub fn is_polynomial(&self, tolerance: T) -> bool {
        // Check if denominator is approximately constant
        let coeffs = self.denominator.unscaled_coefficients();
        if coeffs.is_empty() {
            return false;
        }

        let first = coeffs[0];
        coeffs.iter().all(|&c| (c - first).abs() < tolerance)
    }

    /// Returns the degree of the rational function as an integer.
    ///
    /// Degree of R = degree(P) - degree(Q)
    pub fn rational_degree(&self) -> isize {
        self.numerator.degree() as isize - self.denominator.degree() as isize
    }
}

impl<T: Float + FromPrimitive + fmt::Debug + fmt::Display + std::iter::Sum> fmt::Display for RationalFunction<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "RationalFunction(degrees=[{}, {}], interval=[{}, {}])",
            self.numerator.degree(),
            self.denominator.degree(),
            self.interval().0,
            self.interval().1
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_create_rational_function() {
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 2.0], 0.0, 1.0).unwrap();
        let q = BernsteinPolynomial::from_unscaled(vec![1.0, 1.0], 0.0, 1.0).unwrap();
        let r = RationalFunction::new(p, q).unwrap();

        assert_eq!(r.degrees(), (1, 1));
        assert_eq!(r.interval(), (0.0, 1.0));
    }

    #[test]
    fn test_from_polynomial() {
        // R(x) = x on [0,1]
        let p = BernsteinPolynomial::from_unscaled(vec![0.0, 1.0], 0.0, 1.0).unwrap();
        let r = RationalFunction::from_polynomial(p).unwrap();

        assert_eq!(r.degrees(), (1, 0));
        assert_relative_eq!(r.evaluate(0.5).unwrap(), 0.5, epsilon = 1e-10);
    }

    #[test]
    fn test_evaluate_simple() {
        // R(x) = (1 + 2x) / (1 + x)
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 3.0], 0.0, 1.0).unwrap();
        let q = BernsteinPolynomial::from_unscaled(vec![1.0, 2.0], 0.0, 1.0).unwrap();
        let r = RationalFunction::new(p, q).unwrap();

        // At x=0: (1+0)/(1+0) = 1
        assert_relative_eq!(r.evaluate(0.0).unwrap(), 1.0, epsilon = 1e-10);

        // At x=1: (1+2)/(1+1) = 3/2 = 1.5
        assert_relative_eq!(r.evaluate(1.0).unwrap(), 1.5, epsilon = 1e-10);
    }

    #[test]
    fn test_evaluate_at_pole() {
        // R(x) = x / (x - 0.5) has pole at x = 0.5
        let p = BernsteinPolynomial::from_unscaled(vec![0.0, 1.0], 0.0, 1.0).unwrap();
        let q = BernsteinPolynomial::from_unscaled(vec![0.5, -0.5], 0.0, 1.0).unwrap();
        let r = RationalFunction::new(p, q).unwrap();

        // Should return error at pole
        assert!(r.evaluate(0.5).is_err());
        assert!(r.is_near_pole(0.5, 1e-5));
    }

    #[test]
    fn test_derivative_polynomial() {
        // R(x) = x², R'(x) = 2x
        let p = BernsteinPolynomial::from_unscaled(vec![0.0, 0.0, 1.0], 0.0, 1.0).unwrap();
        let r = RationalFunction::from_polynomial(p).unwrap();
        let dr = r.derivative().unwrap();

        assert_relative_eq!(dr.evaluate(0.0).unwrap(), 0.0, epsilon = 1e-10);
        assert_relative_eq!(dr.evaluate(0.5).unwrap(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(dr.evaluate(1.0).unwrap(), 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_derivative_rational() {
        // R(x) = 1/x on [0.1, 2]
        // R'(x) = -1/x²
        let p = BernsteinPolynomial::constant(1.0, 0.1, 2.0).unwrap();
        let q = BernsteinPolynomial::from_unscaled(vec![0.1, 2.0], 0.1, 2.0).unwrap();
        let r = RationalFunction::new(p, q).unwrap();
        let dr = r.derivative().unwrap();

        // At x=1: R'(1) = -1
        assert_relative_eq!(dr.evaluate(1.0).unwrap(), -1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_constant() {
        let r = RationalFunction::constant(5.0, 0.0, 1.0).unwrap();
        assert_relative_eq!(r.evaluate(0.0).unwrap(), 5.0, epsilon = 1e-10);
        assert_relative_eq!(r.evaluate(0.5).unwrap(), 5.0, epsilon = 1e-10);
        assert_relative_eq!(r.evaluate(1.0).unwrap(), 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_is_polynomial() {
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 2.0, 3.0], 0.0, 1.0).unwrap();
        let r = RationalFunction::from_polynomial(p).unwrap();
        assert!(r.is_polynomial(1e-10));
    }

    #[test]
    fn test_rational_degree() {
        // Degree 2 / Degree 1 = rational degree 1
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 2.0, 3.0], 0.0, 1.0).unwrap();
        let q = BernsteinPolynomial::from_unscaled(vec![1.0, 1.0], 0.0, 1.0).unwrap();
        let r = RationalFunction::new(p, q).unwrap();
        assert_eq!(r.rational_degree(), 1);
    }

    #[test]
    fn test_different_intervals_error() {
        let p = BernsteinPolynomial::constant(1.0, 0.0, 1.0).unwrap();
        let q = BernsteinPolynomial::constant(1.0, 0.0, 2.0).unwrap();
        assert!(RationalFunction::new(p, q).is_err());
    }

    #[test]
    fn test_zero_denominator_error() {
        let p = BernsteinPolynomial::constant(1.0, 0.0, 1.0).unwrap();
        let q = BernsteinPolynomial::constant(0.0, 0.0, 1.0).unwrap();
        assert!(RationalFunction::new(p, q).is_err());
    }
}
