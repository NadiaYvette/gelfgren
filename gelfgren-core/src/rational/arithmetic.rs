//! Arithmetic operations on rational functions.

use super::function::RationalFunction;
use crate::error::Result;
use num_traits::{Float, FromPrimitive};
use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

impl<T: Float + FromPrimitive + fmt::Debug + std::iter::Sum> Add for RationalFunction<T> {
    type Output = Result<Self>;

    /// Adds two rational functions.
    ///
    /// P₁/Q₁ + P₂/Q₂ = (P₁Q₂ + P₂Q₁) / (Q₁Q₂)
    fn add(self, rhs: Self) -> Self::Output {
        if self.interval() != rhs.interval() {
            return Err(crate::error::GelfgrenError::InvalidArgument(
                "Cannot add rational functions on different intervals".to_string(),
            ));
        }

        // (P₁Q₂ + P₂Q₁)
        let term1 = (self.numerator().clone() * rhs.denominator().clone())?;
        let term2 = (rhs.numerator().clone() * self.denominator().clone())?;
        let new_numerator = (term1 + term2)?;

        // Q₁Q₂
        let new_denominator = (self.denominator().clone() * rhs.denominator().clone())?;

        RationalFunction::new(new_numerator, new_denominator)
    }
}

impl<T: Float + FromPrimitive + fmt::Debug + std::iter::Sum> Sub for RationalFunction<T> {
    type Output = Result<Self>;

    /// Subtracts two rational functions.
    ///
    /// P₁/Q₁ - P₂/Q₂ = (P₁Q₂ - P₂Q₁) / (Q₁Q₂)
    fn sub(self, rhs: Self) -> Self::Output {
        if self.interval() != rhs.interval() {
            return Err(crate::error::GelfgrenError::InvalidArgument(
                "Cannot subtract rational functions on different intervals".to_string(),
            ));
        }

        let term1 = (self.numerator().clone() * rhs.denominator().clone())?;
        let term2 = (rhs.numerator().clone() * self.denominator().clone())?;
        let new_numerator = (term1 - term2)?;

        let new_denominator = (self.denominator().clone() * rhs.denominator().clone())?;

        RationalFunction::new(new_numerator, new_denominator)
    }
}

impl<T: Float + FromPrimitive + fmt::Debug + std::iter::Sum> Mul for RationalFunction<T> {
    type Output = Result<Self>;

    /// Multiplies two rational functions.
    ///
    /// (P₁/Q₁) × (P₂/Q₂) = (P₁P₂) / (Q₁Q₂)
    fn mul(self, rhs: Self) -> Self::Output {
        if self.interval() != rhs.interval() {
            return Err(crate::error::GelfgrenError::InvalidArgument(
                "Cannot multiply rational functions on different intervals".to_string(),
            ));
        }

        let new_numerator = (self.numerator().clone() * rhs.numerator().clone())?;
        let new_denominator = (self.denominator().clone() * rhs.denominator().clone())?;

        RationalFunction::new(new_numerator, new_denominator)
    }
}

impl<T: Float + FromPrimitive + fmt::Debug + std::iter::Sum> Neg for RationalFunction<T> {
    type Output = Self;

    /// Negates a rational function.
    ///
    /// -(P/Q) = (-P)/Q
    fn neg(self) -> Self::Output {
        let new_numerator = -self.numerator().clone();
        RationalFunction::new(new_numerator, self.denominator().clone())
            .expect("Negation preserves validity")
    }
}

impl<T: Float + FromPrimitive + fmt::Debug + std::iter::Sum> RationalFunction<T> {
    /// Divides two rational functions.
    ///
    /// (P₁/Q₁) ÷ (P₂/Q₂) = (P₁Q₂) / (Q₁P₂)
    pub fn divide(&self, rhs: &Self) -> Result<Self> {
        if self.interval() != rhs.interval() {
            return Err(crate::error::GelfgrenError::InvalidArgument(
                "Cannot divide rational functions on different intervals".to_string(),
            ));
        }

        let new_numerator = (self.numerator().clone() * rhs.denominator().clone())?;
        let new_denominator = (self.denominator().clone() * rhs.numerator().clone())?;

        Self::new(new_numerator, new_denominator)
    }

    /// Multiplies the rational function by a scalar.
    pub fn scale(&self, scalar: T) -> Self {
        let new_numerator = self.numerator().scale(scalar);
        Self::new(new_numerator, self.denominator().clone())
            .expect("Scaling preserves validity")
    }

    /// Raises the rational function to an integer power.
    ///
    /// R(x)ⁿ = P(x)ⁿ / Q(x)ⁿ
    pub fn pow(&self, n: u32) -> Result<Self> {
        if n == 0 {
            return Self::constant(T::one(), self.interval().0, self.interval().1);
        }
        if n == 1 {
            return Ok(self.clone());
        }

        let mut result = self.clone();
        for _ in 1..n {
            result = (result * self.clone())?;
        }
        Ok(result)
    }

    /// Computes the reciprocal 1/R(x) = Q(x)/P(x).
    pub fn reciprocal(&self) -> Result<Self> {
        Self::new(self.denominator().clone(), self.numerator().clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bernstein::BernsteinPolynomial;
    use approx::assert_relative_eq;

    #[test]
    fn test_add_rational_functions() {
        // 1/x + 1/x = 2/x on [0.1, 2]
        let p = BernsteinPolynomial::constant(1.0, 0.1, 2.0).unwrap();
        let q = BernsteinPolynomial::from_unscaled(vec![0.1, 2.0], 0.1, 2.0).unwrap();
        let r1 = RationalFunction::new(p.clone(), q.clone()).unwrap();
        let r2 = RationalFunction::new(p, q).unwrap();

        let sum = (r1 + r2).unwrap();

        // At x=1: 1/1 + 1/1 = 2
        assert_relative_eq!(sum.evaluate(1.0).unwrap(), 2.0, epsilon = 1e-6);
    }

    #[test]
    fn test_subtract_rational_functions() {
        // 3/x - 1/x = 2/x
        let p1 = BernsteinPolynomial::constant(3.0, 0.1, 2.0).unwrap();
        let p2 = BernsteinPolynomial::constant(1.0, 0.1, 2.0).unwrap();
        let q = BernsteinPolynomial::from_unscaled(vec![0.1, 2.0], 0.1, 2.0).unwrap();

        let r1 = RationalFunction::new(p1, q.clone()).unwrap();
        let r2 = RationalFunction::new(p2, q).unwrap();

        let diff = (r1 - r2).unwrap();

        // At x=1: 3/1 - 1/1 = 2
        assert_relative_eq!(diff.evaluate(1.0).unwrap(), 2.0, epsilon = 1e-6);
    }

    #[test]
    fn test_multiply_rational_functions() {
        // (x/1) × (1/x) = 1
        let p1 = BernsteinPolynomial::from_unscaled(vec![0.1, 2.0], 0.1, 2.0).unwrap();
        let q1 = BernsteinPolynomial::constant(1.0, 0.1, 2.0).unwrap();
        let r1 = RationalFunction::new(p1, q1).unwrap();

        let p2 = BernsteinPolynomial::constant(1.0, 0.1, 2.0).unwrap();
        let q2 = BernsteinPolynomial::from_unscaled(vec![0.1, 2.0], 0.1, 2.0).unwrap();
        let r2 = RationalFunction::new(p2, q2).unwrap();

        let product = (r1 * r2).unwrap();

        // Should be approximately 1 everywhere
        assert_relative_eq!(product.evaluate(0.5).unwrap(), 1.0, epsilon = 1e-4);
        assert_relative_eq!(product.evaluate(1.0).unwrap(), 1.0, epsilon = 1e-4);
        assert_relative_eq!(product.evaluate(1.5).unwrap(), 1.0, epsilon = 1e-4);
    }

    #[test]
    fn test_negate() {
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 2.0], 0.0, 1.0).unwrap();
        let q = BernsteinPolynomial::constant(1.0, 0.0, 1.0).unwrap();
        let r = RationalFunction::new(p, q).unwrap();

        let neg_r = -r.clone();

        assert_relative_eq!(neg_r.evaluate(0.5).unwrap(), -r.evaluate(0.5).unwrap(), epsilon = 1e-10);
    }

    #[test]
    fn test_divide() {
        // (x²/1) ÷ (x/1) = x
        let p1 = BernsteinPolynomial::from_unscaled(vec![0.0, 0.0, 1.0], 0.0, 1.0).unwrap();
        let q1 = BernsteinPolynomial::constant(1.0, 0.0, 1.0).unwrap();
        let r1 = RationalFunction::new(p1, q1).unwrap();

        let p2 = BernsteinPolynomial::from_unscaled(vec![0.0, 1.0], 0.0, 1.0).unwrap();
        let q2 = BernsteinPolynomial::constant(1.0, 0.0, 1.0).unwrap();
        let r2 = RationalFunction::new(p2, q2).unwrap();

        let quotient = r1.divide(&r2).unwrap();

        assert_relative_eq!(quotient.evaluate(0.5).unwrap(), 0.5, epsilon = 1e-6);
        assert_relative_eq!(quotient.evaluate(0.75).unwrap(), 0.75, epsilon = 1e-6);
    }

    #[test]
    fn test_scale() {
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 2.0], 0.0, 1.0).unwrap();
        let q = BernsteinPolynomial::constant(1.0, 0.0, 1.0).unwrap();
        let r = RationalFunction::new(p, q).unwrap();

        let scaled = r.scale(3.0);

        assert_relative_eq!(scaled.evaluate(0.0).unwrap(), 3.0, epsilon = 1e-10);
        assert_relative_eq!(scaled.evaluate(0.5).unwrap(), 4.5, epsilon = 1e-10);
    }

    #[test]
    fn test_pow() {
        // (x)² = x²
        let p = BernsteinPolynomial::from_unscaled(vec![0.0, 1.0], 0.0, 1.0).unwrap();
        let r = RationalFunction::from_polynomial(p).unwrap();
        let r2 = r.pow(2).unwrap();

        assert_relative_eq!(r2.evaluate(0.5).unwrap(), 0.25, epsilon = 1e-6);
        assert_relative_eq!(r2.evaluate(0.75).unwrap(), 0.5625, epsilon = 1e-6);
    }

    #[test]
    fn test_reciprocal() {
        // 1/(x+1) where x+1 is in Bernstein form
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 2.0], 0.0, 1.0).unwrap();
        let r = RationalFunction::from_polynomial(p).unwrap();
        let recip = r.reciprocal().unwrap();

        // At x=0: 1/(0+1) = 1
        assert_relative_eq!(recip.evaluate(0.0).unwrap(), 1.0, epsilon = 1e-10);
        // At x=1: 1/(1+1) = 0.5
        assert_relative_eq!(recip.evaluate(1.0).unwrap(), 0.5, epsilon = 1e-10);
    }
}
