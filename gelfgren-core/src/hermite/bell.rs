//! Bell polynomials for derivative transformations.
//!
//! Implements partial (incomplete) Bell polynomials B_n(x; x_1,...,x_n)
//! which appear in Faà di Bruno's formula and Traub's interpolation formulas.

use num_traits::{Float, FromPrimitive};

/// Represents a Bell polynomial B_n(ω; g_1,...,g_n).
///
/// These polynomials appear in the chain rule for higher derivatives:
/// If h(t) = f(g(t)), then h^(n)(t) involves B_n evaluated at g'(t), g''(t), etc.
pub struct BellPolynomial;

impl BellPolynomial {
    /// Evaluates the partial Bell polynomial B_n(x; x_1, x_2, ..., x_n).
    ///
    /// # Arguments
    ///
    /// * `n` - The order of the Bell polynomial
    /// * `omega` - The main variable ω
    /// * `vars` - Variables [x_1, x_2, ..., x_n] (should have length ≥ n)
    ///
    /// # Formula
    ///
    /// B_n(ω; x_1,...,x_n) = Σ (n!)/(k_1!...k_n!) ∏ᵢ (ωxᵢ/i!)^{kᵢ}
    ///
    /// where sum is over partitions k_1 + 2k_2 + ... + nk_n = n
    ///
    /// # Examples
    ///
    /// - B_0 = 1
    /// - B_1(ω; x_1) = ωx_1
    /// - B_2(ω; x_1, x_2) = ω²x_1² + ωx_2
    /// - B_3(ω; x_1, x_2, x_3) = ω³x_1³ + 3ω²x_1x_2 + ωx_3
    pub fn evaluate<T: Float + FromPrimitive>(n: usize, omega: T, vars: &[T]) -> T {
        if n == 0 {
            return T::one();
        }

        if vars.len() < n {
            return T::zero(); // Not enough variables
        }

        // Use dynamic programming to compute Bell polynomial
        // B_n can be computed recursively: B_n = Σ_{k=1}^n (n-1 choose k-1) x_k B_{n-k}

        let mut bell = vec![T::zero(); n + 1];
        bell[0] = T::one();

        for i in 1..=n {
            let mut sum = T::zero();
            for k in 1..=i {
                let binom = binomial(i - 1, k - 1);
                let binom_f = T::from_usize(binom).unwrap();

                // ω * x_k (no factorial division)
                let term = omega * vars[k - 1];

                sum = sum + binom_f * term * bell[i - k];
            }
            bell[i] = sum;
        }

        bell[n]
    }

    /// Evaluates multiple Bell polynomials B_0, B_1, ..., B_n simultaneously.
    ///
    /// More efficient than calling evaluate() multiple times.
    pub fn evaluate_sequence<T: Float + FromPrimitive>(
        n: usize,
        omega: T,
        vars: &[T],
    ) -> Vec<T> {
        let mut bell = vec![T::zero(); n + 1];
        bell[0] = T::one();

        if n == 0 || vars.is_empty() {
            return bell;
        }

        for i in 1..=n.min(vars.len()) {
            let mut sum = T::zero();
            for k in 1..=i {
                let binom = binomial(i - 1, k - 1);
                let binom_f = T::from_usize(binom).unwrap();

                let term = omega * vars[k - 1];
                sum = sum + binom_f * term * bell[i - k];
            }
            bell[i] = sum;
        }

        bell
    }

    /// Computes the S_r functions used in Traub's interpolation formula.
    ///
    /// S_r(x_i) = (-1)^r (r-1)! Σ_{v≠i} 1/(x_i - x_v)^r
    ///
    /// These appear in the Bell polynomial arguments for Lagrange-Hermite interpolation.
    pub fn compute_s_functions<T: Float + FromPrimitive>(
        points: &[T],
        i: usize,
        max_order: usize,
    ) -> Vec<T> {
        if i >= points.len() {
            return vec![];
        }

        let mut s = vec![T::zero(); max_order + 1];
        let x_i = points[i];

        for r in 1..=max_order {
            let mut sum = T::zero();
            for (v, &x_v) in points.iter().enumerate() {
                if v != i {
                    let diff = x_i - x_v;
                    if diff.abs() > T::from_f64(1e-10).unwrap() {
                        sum = sum + T::one() / diff.powi(r as i32);
                    }
                }
            }

            let sign = if r % 2 == 0 { T::one() } else { -T::one() };
            let factorial_r_minus_1 = T::from_usize(factorial(r - 1)).unwrap();
            s[r] = sign * factorial_r_minus_1 * sum;
        }

        s
    }
}

/// Computes binomial coefficient (n choose k).
fn binomial(n: usize, k: usize) -> usize {
    if k > n {
        return 0;
    }
    if k == 0 || k == n {
        return 1;
    }
    let k = k.min(n - k);
    let mut result = 1;
    for i in 0..k {
        result = result * (n - i) / (i + 1);
    }
    result
}

/// Computes factorial n!
fn factorial(n: usize) -> usize {
    (1..=n).product()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_bell_b0() {
        // B_0 = 1
        let result = BellPolynomial::evaluate(0, 1.0, &[]);
        assert_relative_eq!(result, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_bell_b1() {
        // B_1(ω; x_1) = ωx_1
        let omega = 2.0;
        let vars = vec![3.0];
        let result = BellPolynomial::evaluate(1, omega, &vars);
        assert_relative_eq!(result, 6.0, epsilon = 1e-10); // 2 * 3 = 6
    }

    #[test]
    fn test_bell_b2() {
        // B_2(ω; x_1, x_2) = ω²x_1² + ωx_2
        let omega = 1.0;
        let vars = vec![2.0, 3.0];
        let result = BellPolynomial::evaluate(2, omega, &vars);
        // 1² * 2² + 1 * 3 = 4 + 3 = 7
        assert_relative_eq!(result, 7.0, epsilon = 1e-10);
    }

    #[test]
    fn test_bell_b3() {
        // B_3(ω; x_1, x_2, x_3) = ω³x_1³ + 3ω²x_1x_2 + ωx_3
        let omega = 1.0;
        let vars = vec![1.0, 2.0, 3.0];
        let result = BellPolynomial::evaluate(3, omega, &vars);
        // 1³*1³ + 3*1²*1*2 + 1*3 = 1 + 6 + 3 = 10
        assert_relative_eq!(result, 10.0, epsilon = 1e-10);
    }

    #[test]
    fn test_bell_sequence() {
        let omega = 1.0;
        let vars = vec![1.0, 2.0, 3.0, 4.0];
        let bells = BellPolynomial::evaluate_sequence(3, omega, &vars);

        assert_eq!(bells.len(), 4);
        assert_relative_eq!(bells[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(bells[1], 1.0, epsilon = 1e-10); // ω*x_1 = 1*1 = 1
        assert_relative_eq!(bells[2], 3.0, epsilon = 1e-10); // ω²x_1² + ωx_2 = 1 + 2 = 3
        assert_relative_eq!(bells[3], 10.0, epsilon = 1e-10); // As computed above
    }

    #[test]
    fn test_s_functions() {
        // Points: [0, 1, 2]
        // For i=1 (x_1 = 1):
        // S_1(x_1) = -0! * (1/(1-0) + 1/(1-2)) = -(1 - 1) = 0
        // S_2(x_1) = 1! * (1/(1-0)² + 1/(1-2)²) = 1 + 1 = 2
        let points = vec![0.0, 1.0, 2.0];
        let s = BellPolynomial::compute_s_functions(&points, 1, 2);

        assert_eq!(s.len(), 3);
        assert_relative_eq!(s[1], 0.0, epsilon = 1e-10);
        assert_relative_eq!(s[2], 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_binomial() {
        assert_eq!(binomial(5, 0), 1);
        assert_eq!(binomial(5, 1), 5);
        assert_eq!(binomial(5, 2), 10);
        assert_eq!(binomial(5, 3), 10);
        assert_eq!(binomial(10, 3), 120);
    }

    #[test]
    fn test_factorial() {
        assert_eq!(factorial(0), 1);
        assert_eq!(factorial(1), 1);
        assert_eq!(factorial(5), 120);
        assert_eq!(factorial(10), 3628800);
    }
}
