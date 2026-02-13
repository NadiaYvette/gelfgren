//! Differentiation and integration of Bernstein polynomials.
//!
//! Implements algorithms from Section 5.1 of Farouki & Rajan (1987).

use super::polynomial::BernsteinPolynomial;
use crate::error::Result;
use num_traits::{Float, FromPrimitive};

impl<T: Float + FromPrimitive + std::fmt::Debug + std::iter::Sum> BernsteinPolynomial<T> {
    /// Computes the derivative dP/dx.
    ///
    /// Uses equation (57) from Farouki-Rajan:
    /// D̂_k^{n-1} = n/(b-a) [Ĉ_{k+1}^n / C(n,k+1) - Ĉ_k^n / C(n,k)] C(n-1,k)
    ///
    /// Simplified: d/dx P(x) has degree n-1 with scaled coefficients
    /// D̂_k = n/(b-a) * (C_{k+1} - C_k) * C(n-1,k)
    ///
    /// where C_k are the unscaled coefficients.
    pub fn derivative(&self) -> Result<Self> {
        let n = self.degree();
        let (a, b) = self.interval();

        if n == 0 {
            // Derivative of constant is zero
            return Self::zero_on_interval(a, b);
        }
        let interval_length = b - a;

        // Convert to unscaled coefficients
        let unscaled: Vec<T> = self.unscaled_coefficients();

        // Compute unscaled derivative coefficients
        let mut d_unscaled = vec![T::zero(); n];
        let scale = T::from_usize(n).unwrap() / interval_length;

        for k in 0..n {
            d_unscaled[k] = scale * (unscaled[k + 1] - unscaled[k]);
        }

        // Convert back to scaled form
        Self::from_unscaled(d_unscaled, a, b)
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

    /// Computes the indefinite integral ∫P(x)dx.
    ///
    /// Uses equation (59) from Farouki-Rajan:
    /// Î_k^{n+1} = (b-a)/(n+1) Σ_{j=0}^{k-1} Ĉ_j^n / C(n,j) * C(n+1,k)
    ///
    /// The constant of integration is chosen so that ∫_a^x P(t)dt = 0 at x=a.
    pub fn integral(&self) -> Result<Self> {
        let n = self.degree();
        let (a, b) = self.interval();
        let interval_length = b - a;

        // Convert to unscaled coefficients
        let unscaled: Vec<T> = self.unscaled_coefficients();

        // Compute unscaled integral coefficients
        let mut i_unscaled = vec![T::zero(); n + 2];
        let scale = interval_length / T::from_usize(n + 1).unwrap();

        // First coefficient is 0 (integral at left endpoint)
        i_unscaled[0] = T::zero();

        // Cumulative sum for integral
        for k in 1..=n + 1 {
            i_unscaled[k] = i_unscaled[k - 1] + scale * unscaled[k - 1];
        }

        Self::from_unscaled(i_unscaled, a, b)
    }

    /// Computes the definite integral ∫_a^b P(x)dx.
    ///
    /// Uses the fact that for Bernstein polynomials on [a,b]:
    /// ∫_a^b P(x)dx = (b-a) Σ_{k=0}^n C_k^n / (n+1)
    /// where C_k^n are unscaled coefficients.
    pub fn definite_integral(&self) -> T {
        let n = self.degree();
        let (a, b) = self.interval();
        let interval_length = b - a;

        let unscaled = self.unscaled_coefficients();
        let sum: T = unscaled.iter().copied().sum();

        interval_length * sum / T::from_usize(n + 1).unwrap()
    }

    /// Computes the definite integral over a subinterval [x1, x2] ⊆ [a, b].
    pub fn definite_integral_on(&self, x1: T, x2: T) -> T {
        let indefinite = self.integral().expect("Integral computation failed");
        indefinite.evaluate(x2) - indefinite.evaluate(x1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_derivative_constant() {
        let p = BernsteinPolynomial::constant(5.0, 0.0, 1.0).unwrap();
        let dp = p.derivative().unwrap();

        assert_eq!(dp.degree(), 0);
        assert_relative_eq!(dp.evaluate(0.5), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_derivative_linear() {
        // P(x) = 1 + 2x on [0,1]
        // dP/dx = 2
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 3.0], 0.0, 1.0).unwrap();
        let dp = p.derivative().unwrap();

        assert_eq!(dp.degree(), 0);
        assert_relative_eq!(dp.evaluate(0.0), 2.0, epsilon = 1e-10);
        assert_relative_eq!(dp.evaluate(0.5), 2.0, epsilon = 1e-10);
        assert_relative_eq!(dp.evaluate(1.0), 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_derivative_quadratic() {
        // P(x) = x^2 on [0,1]
        // dP/dx = 2x
        let p = BernsteinPolynomial::from_unscaled(vec![0.0, 0.0, 1.0], 0.0, 1.0).unwrap();
        let dp = p.derivative().unwrap();

        assert_eq!(dp.degree(), 1);
        assert_relative_eq!(dp.evaluate(0.0), 0.0, epsilon = 1e-10);
        assert_relative_eq!(dp.evaluate(0.5), 1.0, epsilon = 1e-10);
        assert_relative_eq!(dp.evaluate(1.0), 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_derivative_scaled_interval() {
        // P(x) = x on [0, 2]
        // dP/dx = 1
        let p = BernsteinPolynomial::from_unscaled(vec![0.0, 2.0], 0.0, 2.0).unwrap();
        let dp = p.derivative().unwrap();

        assert_relative_eq!(dp.evaluate(0.0), 1.0, epsilon = 1e-10);
        assert_relative_eq!(dp.evaluate(1.0), 1.0, epsilon = 1e-10);
        assert_relative_eq!(dp.evaluate(2.0), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_second_derivative() {
        // P(x) = x^2
        // dP/dx = 2x
        // d²P/dx² = 2
        let p = BernsteinPolynomial::from_unscaled(vec![0.0, 0.0, 1.0], 0.0, 1.0).unwrap();
        let d2p = p.derivative_n(2).unwrap();

        assert_eq!(d2p.degree(), 0);
        assert_relative_eq!(d2p.evaluate(0.5), 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_integral_constant() {
        // P(x) = 5
        // ∫P(x)dx = 5x (with constant 0)
        let p = BernsteinPolynomial::constant(5.0, 0.0, 1.0).unwrap();
        let ip = p.integral().unwrap();

        assert_eq!(ip.degree(), 1);
        assert_relative_eq!(ip.evaluate(0.0), 0.0, epsilon = 1e-10);
        assert_relative_eq!(ip.evaluate(0.5), 2.5, epsilon = 1e-10);
        assert_relative_eq!(ip.evaluate(1.0), 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_integral_linear() {
        // P(x) = 2x on [0,1]
        // ∫P(x)dx = x^2
        let p = BernsteinPolynomial::from_unscaled(vec![0.0, 2.0], 0.0, 1.0).unwrap();
        let ip = p.integral().unwrap();

        assert_relative_eq!(ip.evaluate(0.0), 0.0, epsilon = 1e-10);
        assert_relative_eq!(ip.evaluate(0.5), 0.25, epsilon = 1e-10);
        assert_relative_eq!(ip.evaluate(1.0), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_definite_integral_constant() {
        // ∫_0^1 5 dx = 5
        let p = BernsteinPolynomial::constant(5.0, 0.0, 1.0).unwrap();
        assert_relative_eq!(p.definite_integral(), 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_definite_integral_linear() {
        // ∫_0^1 (1 + 2x) dx = [x + x^2]_0^1 = 2
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 3.0], 0.0, 1.0).unwrap();
        assert_relative_eq!(p.definite_integral(), 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_definite_integral_quadratic() {
        // ∫_0^1 x^2 dx = 1/3
        let p = BernsteinPolynomial::from_unscaled(vec![0.0, 0.0, 1.0], 0.0, 1.0).unwrap();
        assert_relative_eq!(p.definite_integral(), 1.0 / 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_fundamental_theorem() {
        // Verify: d/dx ∫P(x)dx = P(x)
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 2.0, 3.0], 0.0, 1.0).unwrap();
        let ip = p.integral().unwrap();
        let dip = ip.derivative().unwrap();

        for x in [0.0, 0.25, 0.5, 0.75, 1.0] {
            assert_relative_eq!(p.evaluate(x), dip.evaluate(x), epsilon = 1e-9);
        }
    }

    #[test]
    fn test_definite_integral_on_subinterval() {
        // P(x) = x on [0,2]
        // ∫_0.5^1.5 x dx = [x²/2]_0.5^1.5 = 1.125 - 0.125 = 1.0
        let p = BernsteinPolynomial::from_unscaled(vec![0.0, 2.0], 0.0, 2.0).unwrap();
        assert_relative_eq!(p.definite_integral_on(0.5, 1.5), 1.0, epsilon = 1e-10);
    }
}
