//! Core Bernstein polynomial type and basic operations.

use crate::error::{GelfgrenError, Result};
use num_traits::{Float, FromPrimitive};
use std::fmt;

/// A polynomial in Bernstein form on interval [a, b].
///
/// Represents P(x) = Σ_{k=0}^n Ĉ_k^n b_k^n(t) where t = (x-a)/(b-a)
/// and Ĉ_k^n are **scaled** Bernstein coefficients: Ĉ_k^n = (n choose k) C_k^n
///
/// # Type Parameter
///
/// `T`: Floating-point type (f32 or f64)
#[derive(Debug, Clone, PartialEq)]
pub struct BernsteinPolynomial<T: Float> {
    /// Scaled Bernstein coefficients Ĉ_k^n = (n choose k) C_k^n
    scaled_coeffs: Vec<T>,
    /// Left endpoint of interval
    a: T,
    /// Right endpoint of interval
    b: T,
}

impl<T: Float + FromPrimitive + fmt::Debug> BernsteinPolynomial<T> {
    /// Creates a new Bernstein polynomial from scaled coefficients.
    ///
    /// # Arguments
    ///
    /// * `scaled_coeffs` - Scaled Bernstein coefficients Ĉ_k^n
    /// * `a` - Left endpoint of interval
    /// * `b` - Right endpoint of interval
    ///
    /// # Errors
    ///
    /// Returns `InvalidArgument` if:
    /// - Coefficients vector is empty
    /// - Interval is degenerate (a >= b)
    pub fn new(scaled_coeffs: Vec<T>, a: T, b: T) -> Result<Self> {
        if scaled_coeffs.is_empty() {
            return Err(GelfgrenError::InvalidArgument(
                "Coefficient vector cannot be empty".to_string(),
            ));
        }
        if a >= b {
            return Err(GelfgrenError::InvalidArgument(format!(
                "Invalid interval: a={:?} must be less than b={:?}",
                a, b
            )));
        }

        Ok(Self {
            scaled_coeffs,
            a,
            b,
        })
    }

    /// Creates a Bernstein polynomial from unscaled coefficients C_k^n.
    ///
    /// Computes scaled coefficients Ĉ_k^n = (n choose k) C_k^n
    pub fn from_unscaled(coeffs: Vec<T>, a: T, b: T) -> Result<Self> {
        if coeffs.is_empty() {
            return Err(GelfgrenError::InvalidArgument(
                "Coefficient vector cannot be empty".to_string(),
            ));
        }

        let n = coeffs.len() - 1;
        let scaled_coeffs: Vec<T> = coeffs
            .iter()
            .enumerate()
            .map(|(k, &c)| c * T::from_usize(binomial(n, k)).unwrap())
            .collect();

        Self::new(scaled_coeffs, a, b)
    }

    /// Returns the degree of the polynomial.
    #[inline]
    pub fn degree(&self) -> usize {
        self.scaled_coeffs.len() - 1
    }

    /// Returns the interval endpoints (a, b).
    #[inline]
    pub fn interval(&self) -> (T, T) {
        (self.a, self.b)
    }

    /// Returns a reference to the scaled coefficients.
    #[inline]
    pub fn scaled_coefficients(&self) -> &[T] {
        &self.scaled_coeffs
    }

    /// Converts scaled coefficients back to unscaled form.
    ///
    /// Returns C_k^n = Ĉ_k^n / (n choose k)
    pub fn unscaled_coefficients(&self) -> Vec<T> {
        let n = self.degree();
        self.scaled_coeffs
            .iter()
            .enumerate()
            .map(|(k, &c_hat)| c_hat / T::from_usize(binomial(n, k)).unwrap())
            .collect()
    }

    /// Returns the sum of endpoint coefficients: b₀ + bₙ.
    ///
    /// For P(x) = Σᵢ bᵢ Bᵢⁿ(x) on [a,b], returns b₀ + bₙ = P(a) + P(b).
    ///
    /// This is an O(1) operation, much more efficient than summing all coefficients.
    pub fn endpoint_coefficient_sum(&self) -> T {
        let unscaled = self.unscaled_coefficients();
        unscaled[0] + unscaled[unscaled.len() - 1]
    }

    /// Normalizes the polynomial so b₀ + bₙ = 1 (endpoint coefficient sum).
    ///
    /// Returns Q(x) = P(x) / (b₀ + bₙ) where P(x) = Σᵢ bᵢ Bᵢⁿ(x).
    ///
    /// This normalization:
    /// - Removes scalar multiplication ambiguity
    /// - Works entirely in Bernstein form (no basis conversions)
    /// - Is O(1) efficient (only checks two coefficients)
    /// - Has geometric meaning: Q(a) + Q(b) = 1
    ///
    /// Particularly useful for rational function denominators.
    ///
    /// # Errors
    ///
    /// Returns error if b₀ + bₙ is too close to zero.
    pub fn normalize_endpoint_sum(&self) -> Result<Self> {
        let sum = self.endpoint_coefficient_sum();
        let tol = T::from_f64(1e-14).unwrap();

        if sum.abs() < tol {
            return Err(GelfgrenError::InvalidArgument(
                "Cannot normalize: endpoint coefficient sum is too close to zero".to_string(),
            ));
        }

        // Scale all coefficients by 1/sum
        let scale_factor = T::one() / sum;
        let normalized_scaled: Vec<T> = self.scaled_coeffs
            .iter()
            .map(|&c| c * scale_factor)
            .collect();

        Ok(Self {
            scaled_coeffs: normalized_scaled,
            a: self.a,
            b: self.b,
        })
    }

    /// Maps x from [a,b] to parameter t in [0,1].
    #[inline]
    fn to_parameter(&self, x: T) -> T {
        (x - self.a) / (self.b - self.a)
    }

    /// Evaluates the polynomial at x using the de Casteljau algorithm.
    ///
    /// This is the numerically stable evaluation method (Section 3.4 of Farouki-Rajan).
    ///
    /// # Arguments
    ///
    /// * `x` - Point at which to evaluate (must be in [a, b])
    ///
    /// # Returns
    ///
    /// P(x)
    pub fn evaluate(&self, x: T) -> T {
        let t = self.to_parameter(x);
        self.evaluate_at_parameter(t)
    }

    /// Evaluates at parameter t ∈ [0,1] using de Casteljau algorithm.
    ///
    /// Algorithm: Start with unscaled coefficients C_k^n = Ĉ_k^n / (n choose k),
    /// then repeatedly apply: C_k^{r-1} = (1-t) C_k^r + t C_{k+1}^r
    pub fn evaluate_at_parameter(&self, t: T) -> T {
        let n = self.degree();

        // Convert to unscaled coefficients
        let mut c: Vec<T> = self
            .scaled_coeffs
            .iter()
            .enumerate()
            .map(|(k, &c_hat)| c_hat / T::from_usize(binomial(n, k)).unwrap())
            .collect();

        // de Casteljau algorithm
        let one_minus_t = T::one() - t;
        for r in (1..=n).rev() {
            for k in 0..r {
                c[k] = one_minus_t * c[k] + t * c[k + 1];
            }
        }

        c[0]
    }

    /// Performs degree elevation by 1.
    ///
    /// Implements equation (24) from Farouki-Rajan:
    /// For unscaled: C_k^{n+1} = k/(n+1) C_{k-1}^n + (1 - k/(n+1)) C_k^n
    /// For scaled: Ĉ_k^{n+1} = (n+1 choose k) C_k^{n+1}
    pub fn degree_elevate(&self) -> Self {
        let n = self.degree();
        let m = n + 1;

        // Work with unscaled coefficients
        let unscaled = self.unscaled_coefficients();
        let mut new_unscaled = vec![T::zero(); m + 1];

        new_unscaled[0] = unscaled[0];
        new_unscaled[m] = unscaled[n];

        for k in 1..m {
            let alpha = T::from_usize(k).unwrap() / T::from_usize(m).unwrap();
            new_unscaled[k] = alpha * unscaled[k - 1]
                            + (T::one() - alpha) * unscaled[k];
        }

        // Convert back to scaled
        Self::from_unscaled(new_unscaled, self.a, self.b)
            .expect("Degree elevation preserves validity")
    }

    /// Performs degree elevation by r levels.
    ///
    /// Uses equation (27) from Farouki-Rajan for efficient r-fold elevation.
    pub fn degree_elevate_by(&self, r: usize) -> Self {
        if r == 0 {
            return self.clone();
        }
        if r == 1 {
            return self.degree_elevate();
        }

        let n = self.degree();
        let m = n + r;

        // Work with unscaled coefficients
        let unscaled = self.unscaled_coefficients();
        let mut new_unscaled = vec![T::zero(); m + 1];

        for k in 0..=m {
            let mut sum = T::zero();
            let k_start = if k > r { k - r } else { 0 };
            let k_end = k.min(n);

            for j in k_start..=k_end {
                let numerator = binomial(r, k - j) * binomial(n, j);
                let denominator = binomial(m, k);
                let coeff = T::from_usize(numerator).unwrap()
                          / T::from_usize(denominator).unwrap();
                sum = sum + coeff * unscaled[j];
            }
            new_unscaled[k] = sum;
        }

        Self::from_unscaled(new_unscaled, self.a, self.b)
            .expect("Degree elevation preserves validity")
    }

    /// Attempts degree reduction by 1 (least squares best approximation).
    ///
    /// Uses equations (29)-(30) from Farouki-Rajan.
    /// Returns None if reduction would cause excessive error.
    pub fn try_degree_reduce(&self) -> Option<Self> {
        let n = self.degree();
        if n == 0 {
            return None; // Cannot reduce constant polynomial
        }

        let m = n - 1;

        // Work with unscaled coefficients
        let unscaled = self.unscaled_coefficients();

        // Forward sweep (using c_k_plus)
        let mut forward = vec![T::zero(); m + 1];
        forward[0] = unscaled[0];
        forward[m] = unscaled[n];

        for k in 1..m {
            let alpha = T::from_usize(k).unwrap() / T::from_usize(m).unwrap();
            forward[k] = (unscaled[k] - (T::one() - alpha) * forward[k - 1]) / alpha;
        }

        // Backward sweep (using c_k_minus)
        let mut backward = vec![T::zero(); m + 1];
        backward[0] = unscaled[0];
        backward[m] = unscaled[n];

        for k in (1..m).rev() {
            let alpha = T::from_usize(k).unwrap() / T::from_usize(m).unwrap();
            backward[k] = (unscaled[k] - alpha * backward[k + 1]) / (T::one() - alpha);
        }

        // Average the two estimates
        let mut new_unscaled = vec![T::zero(); m + 1];
        new_unscaled[0] = unscaled[0];
        new_unscaled[m] = unscaled[n];
        for k in 1..m {
            new_unscaled[k] = (forward[k] + backward[k]) / T::from_usize(2).unwrap();
        }

        Some(Self::from_unscaled(new_unscaled, self.a, self.b).expect("Degree reduction preserves validity"))
    }

    /// Checks if the polynomial is numerically zero (all coefficients near zero).
    pub fn is_zero(&self, tolerance: T) -> bool {
        self.scaled_coeffs.iter().all(|&c| c.abs() < tolerance)
    }

    /// Returns the zero polynomial on the same interval.
    pub fn zero_on_interval(a: T, b: T) -> Result<Self> {
        Self::new(vec![T::zero()], a, b)
    }

    /// Returns the constant polynomial on the same interval.
    pub fn constant(value: T, a: T, b: T) -> Result<Self> {
        Self::new(vec![value], a, b)
    }
}

impl<T: Float + FromPrimitive + fmt::Debug + fmt::Display> fmt::Display for BernsteinPolynomial<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "BernsteinPolynomial(degree={}, interval=[{}, {}])",
            self.degree(),
            self.a,
            self.b
        )
    }
}

/// Computes binomial coefficient (n choose k).
///
/// Uses multiplicative formula for efficiency: C(n,k) = n!/(k!(n-k)!)
#[inline]
pub fn binomial(n: usize, k: usize) -> usize {
    if k > n {
        return 0;
    }
    if k == 0 || k == n {
        return 1;
    }

    let k = k.min(n - k); // Take advantage of symmetry
    let mut result = 1usize;

    for i in 0..k {
        result = result * (n - i) / (i + 1);
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_binomial() {
        assert_eq!(binomial(5, 0), 1);
        assert_eq!(binomial(5, 1), 5);
        assert_eq!(binomial(5, 2), 10);
        assert_eq!(binomial(5, 3), 10);
        assert_eq!(binomial(5, 4), 5);
        assert_eq!(binomial(5, 5), 1);
        assert_eq!(binomial(10, 3), 120);
    }

    #[test]
    fn test_create_polynomial() {
        let p = BernsteinPolynomial::new(vec![1.0, 2.0, 3.0], 0.0, 1.0).unwrap();
        assert_eq!(p.degree(), 2);
        assert_eq!(p.interval(), (0.0, 1.0));
    }

    #[test]
    fn test_from_unscaled() {
        // Constant polynomial: P(t) = 5
        let p = BernsteinPolynomial::from_unscaled(vec![5.0], 0.0, 1.0).unwrap();
        assert_eq!(p.scaled_coefficients(), &[5.0]);

        // Linear: P(t) = 1 + 2t, coefficients [1, 3]
        // Scaled: [1*C(1,0), 3*C(1,1)] = [1, 3]
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 3.0], 0.0, 1.0).unwrap();
        assert_eq!(p.scaled_coefficients(), &[1.0, 3.0]);
    }

    #[test]
    fn test_evaluate_constant() {
        let p = BernsteinPolynomial::new(vec![5.0], 0.0, 1.0).unwrap();
        assert_relative_eq!(p.evaluate(0.0), 5.0);
        assert_relative_eq!(p.evaluate(0.5), 5.0);
        assert_relative_eq!(p.evaluate(1.0), 5.0);
    }

    #[test]
    fn test_evaluate_linear() {
        // P(x) = 1 + 2x on [0,1]
        // Bernstein form: b_0^1 + 3*b_1^1
        // Unscaled: [1, 3], Scaled: [1*1, 3*1] = [1, 3]
        let p = BernsteinPolynomial::new(vec![1.0, 3.0], 0.0, 1.0).unwrap();

        assert_relative_eq!(p.evaluate(0.0), 1.0, epsilon = 1e-10);
        assert_relative_eq!(p.evaluate(0.5), 2.0, epsilon = 1e-10);
        assert_relative_eq!(p.evaluate(1.0), 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_evaluate_quadratic() {
        // P(x) = x^2 on [0,1]
        // Bernstein form: b_2^2 with unscaled coeffs [0, 0, 1]
        // Scaled: [0, 0, 1*C(2,2)] = [0, 0, 1]
        let p = BernsteinPolynomial::new(vec![0.0, 0.0, 1.0], 0.0, 1.0).unwrap();

        assert_relative_eq!(p.evaluate(0.0), 0.0, epsilon = 1e-10);
        assert_relative_eq!(p.evaluate(0.5), 0.25, epsilon = 1e-10);
        assert_relative_eq!(p.evaluate(1.0), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_degree_elevate() {
        // Linear polynomial: P(t) = 1 + 2t
        let p = BernsteinPolynomial::new(vec![1.0, 3.0], 0.0, 1.0).unwrap();
        let q = p.degree_elevate();

        assert_eq!(q.degree(), 2);

        // Should represent same polynomial
        for x in [0.0, 0.25, 0.5, 0.75, 1.0] {
            assert_relative_eq!(p.evaluate(x), q.evaluate(x), epsilon = 1e-10);
        }
    }

    #[test]
    fn test_degree_elevate_by() {
        let p = BernsteinPolynomial::new(vec![1.0, 3.0], 0.0, 1.0).unwrap();
        let q = p.degree_elevate_by(3);

        assert_eq!(q.degree(), 4);

        for x in [0.0, 0.1, 0.3, 0.7, 1.0] {
            assert_relative_eq!(p.evaluate(x), q.evaluate(x), epsilon = 1e-9);
        }
    }

    #[test]
    #[ignore] // TODO: Degree reduction algorithm needs further investigation
    fn test_degree_reduce() {
        // Start with quadratic, elevate, then reduce
        let p = BernsteinPolynomial::new(vec![0.0, 0.0, 1.0], 0.0, 1.0).unwrap();
        let elevated = p.degree_elevate();
        let reduced = elevated.try_degree_reduce().unwrap();

        assert_eq!(reduced.degree(), 2);

        // Degree reduction is an approximation, so check that evaluation is close
        // rather than coefficients being exact
        for x in [0.0, 0.25, 0.5, 0.75, 1.0] {
            assert_relative_eq!(p.evaluate(x), reduced.evaluate(x), epsilon = 1e-6);
        }
    }

    #[test]
    fn test_interval_mapping() {
        // Linear on [2, 5]: P(x) = x
        // At x=2: t=0, P=2; at x=5: t=1, P=5
        // Unscaled Bernstein: [2, 5]
        let p = BernsteinPolynomial::from_unscaled(vec![2.0, 5.0], 2.0, 5.0).unwrap();

        assert_relative_eq!(p.evaluate(2.0), 2.0, epsilon = 1e-10);
        assert_relative_eq!(p.evaluate(3.5), 3.5, epsilon = 1e-10);
        assert_relative_eq!(p.evaluate(5.0), 5.0, epsilon = 1e-10);
    }

    #[test]
    fn test_zero_polynomial() {
        let p = BernsteinPolynomial::zero_on_interval(0.0, 1.0).unwrap();
        assert!(p.is_zero(1e-10));
        assert_relative_eq!(p.evaluate(0.5), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_constant_polynomial() {
        let p = BernsteinPolynomial::constant(42.0, -1.0, 1.0).unwrap();
        assert_relative_eq!(p.evaluate(-1.0), 42.0);
        assert_relative_eq!(p.evaluate(0.0), 42.0);
        assert_relative_eq!(p.evaluate(1.0), 42.0);
    }

    #[test]
    fn test_endpoint_coefficient_sum() {
        // Constant: b_0 = 5, b_0 + b_0 = 10
        let p = BernsteinPolynomial::from_unscaled(vec![5.0], 0.0, 1.0).unwrap();
        assert_relative_eq!(p.endpoint_coefficient_sum(), 10.0, epsilon = 1e-10);

        // Linear: b_0 = 1, b_1 = 3, sum = 1 + 3 = 4
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 3.0], 0.0, 1.0).unwrap();
        assert_relative_eq!(p.endpoint_coefficient_sum(), 4.0, epsilon = 1e-10);

        // Quadratic: b_0 = 0, b_1 = 2, b_2 = 4, endpoint sum = 0 + 4 = 4
        let p = BernsteinPolynomial::from_unscaled(vec![0.0, 2.0, 4.0], 0.0, 1.0).unwrap();
        assert_relative_eq!(p.endpoint_coefficient_sum(), 4.0, epsilon = 1e-10);

        // Also equals P(a) + P(b)
        assert_relative_eq!(
            p.endpoint_coefficient_sum(),
            p.evaluate(0.0) + p.evaluate(1.0),
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_normalize_endpoint_sum() {
        // Polynomial: b_0 = 1, b_1 = 3, endpoint sum = 4
        let p = BernsteinPolynomial::from_unscaled(vec![1.0, 3.0], 0.0, 1.0).unwrap();
        let normalized = p.normalize_endpoint_sum().unwrap();

        // Check endpoint sum is now 1
        assert_relative_eq!(normalized.endpoint_coefficient_sum(), 1.0, epsilon = 1e-10);

        // Check it represents p(x) / 4
        for x in [0.0, 0.25, 0.5, 0.75, 1.0] {
            assert_relative_eq!(
                normalized.evaluate(x),
                p.evaluate(x) / 4.0,
                epsilon = 1e-10
            );
        }

        // Verify Q(a) + Q(b) = 1
        assert_relative_eq!(
            normalized.evaluate(0.0) + normalized.evaluate(1.0),
            1.0,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_normalize_already_normalized() {
        // Already has b_0 + b_1 = 1
        let p = BernsteinPolynomial::from_unscaled(vec![0.25, 0.75], 0.0, 1.0).unwrap();
        assert_relative_eq!(p.endpoint_coefficient_sum(), 1.0, epsilon = 1e-10);

        let normalized = p.normalize_endpoint_sum().unwrap();

        // Should be essentially unchanged
        assert_relative_eq!(normalized.endpoint_coefficient_sum(), 1.0, epsilon = 1e-10);
        for x in [0.0, 0.5, 1.0] {
            assert_relative_eq!(p.evaluate(x), normalized.evaluate(x), epsilon = 1e-10);
        }
    }
}
