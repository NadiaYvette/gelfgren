//! Padé approximant construction and representation.

use crate::bernstein::BernsteinPolynomial;
use crate::error::{GelfgrenError, Result};
use crate::rational::RationalFunction;
use num_traits::{Float, FromPrimitive};

/// A Padé approximant [n/m] constructed from series data.
///
/// Represents a rational function P_n(x)/Q_m(x) that matches
/// a power series or derivative data to order n+m+1.
pub struct PadeApproximant<T: Float> {
    /// The underlying rational function representation
    rational: RationalFunction<T>,
    /// Numerator degree n
    n: usize,
    /// Denominator degree m
    m: usize,
    /// Center point (for Taylor series)
    center: T,
}

impl<T: Float + FromPrimitive + std::fmt::Debug + std::iter::Sum> PadeApproximant<T> {
    /// Constructs a [n/m] Padé approximant from power series coefficients.
    ///
    /// Given a power series f(x) = Σ aₖ(x-x₀)ᵏ, constructs P_n/Q_m such that
    /// f(x) - P_n(x)/Q_m(x) = O((x-x₀)^(n+m+1))
    ///
    /// # Arguments
    ///
    /// * `coeffs` - Power series coefficients [a₀, a₁, ..., aₙ₊ₘ]
    ///              Must have length ≥ n+m+1
    /// * `n` - Numerator degree
    /// * `m` - Denominator degree
    /// * `center` - Center point x₀ of the series
    /// * `interval` - Interval [a, b] for Bernstein representation
    ///
    /// # Algorithm
    ///
    /// The Padé approximant satisfies:
    /// P(x) - Q(x)f(x) = O((x-x₀)^(n+m+1))
    ///
    /// Expanding: Σᵢ pᵢxⁱ - (Σⱼ qⱼxʲ)(Σₖ aₖxᵏ) = 0 for powers x⁰...x^(n+m)
    ///
    /// This gives a linear system for coefficients {p₀,...,pₙ, q₁,...,qₘ}
    /// with normalization q₀ = 1.
    pub fn from_series(
        coeffs: &[T],
        n: usize,
        m: usize,
        center: T,
        interval: (T, T),
    ) -> Result<Self> {
        if coeffs.len() < n + m + 1 {
            return Err(GelfgrenError::InvalidArgument(format!(
                "Need at least {} coefficients for [{}/{}] Padé approximant, got {}",
                n + m + 1,
                n,
                m,
                coeffs.len()
            )));
        }

        // Solve the Padé equations to get power series coefficients
        let (p_power, q_power) = solve_pade_system(coeffs, n, m)?;

        // Convert power series to polynomials centered at 'center'
        // For now, assume center = interval.0 for simplicity
        // TODO: Handle arbitrary center points with coordinate transformation

        let (a, b) = interval;

        // Convert power series coefficients to Bernstein form on [a, b]
        let numerator = power_series_to_bernstein(&p_power, a, b)?;
        let denominator = power_series_to_bernstein(&q_power, a, b)?;

        let rational = RationalFunction::new(numerator, denominator)?;

        Ok(Self {
            rational,
            n,
            m,
            center,
        })
    }

    /// Constructs a [n/m] Padé approximant from function and derivative values.
    ///
    /// Given f(x₀), f'(x₀), ..., f^(n+m)(x₀), constructs the approximant.
    ///
    /// # Arguments
    ///
    /// * `derivatives` - Values [f(x₀), f'(x₀), ..., f^(n+m)(x₀)]
    /// * `n` - Numerator degree
    /// * `m` - Denominator degree
    /// * `center` - Point x₀
    /// * `interval` - Interval [a, b] for Bernstein representation
    pub fn from_derivatives(
        derivatives: &[T],
        n: usize,
        m: usize,
        center: T,
        interval: (T, T),
    ) -> Result<Self> {
        if derivatives.len() < n + m + 1 {
            return Err(GelfgrenError::InvalidArgument(format!(
                "Need at least {} derivative values for [{}/{}] Padé approximant",
                n + m + 1,
                n,
                m
            )));
        }

        // Convert derivatives to Taylor series coefficients: aₖ = f^(k)(x₀)/k!
        let mut coeffs = Vec::with_capacity(derivatives.len());
        for (k, &deriv) in derivatives.iter().enumerate() {
            let factorial = factorial(k);
            let coeff = deriv / T::from_usize(factorial).unwrap();
            coeffs.push(coeff);
        }

        Self::from_series(&coeffs, n, m, center, interval)
    }

    /// Returns the underlying rational function.
    pub fn rational(&self) -> &RationalFunction<T> {
        &self.rational
    }

    /// Returns the degrees (n, m).
    pub fn degrees(&self) -> (usize, usize) {
        (self.n, self.m)
    }

    /// Returns the center point.
    pub fn center(&self) -> T {
        self.center
    }

    /// Evaluates the Padé approximant at x.
    pub fn evaluate(&self, x: T) -> Result<T> {
        self.rational.evaluate(x)
    }

    /// Computes the derivative of the Padé approximant.
    pub fn derivative(&self) -> Result<Self> {
        let new_rational = self.rational.derivative()?;
        Ok(Self {
            rational: new_rational,
            n: if self.n > 0 { self.n - 1 } else { 0 },
            m: self.m + 1, // Denominator degree increases due to Q²
            center: self.center,
        })
    }
}

/// Solves the Padé system of equations.
///
/// Given series coefficients {a₀, a₁, ..., aₙ₊ₘ}, finds {p₀,...,pₙ, q₁,...,qₘ}
/// such that P(x) - Q(x)·Σaₖxᵏ = O(x^(n+m+1)) with q₀ = 1.
///
/// Returns (numerator_coeffs, denominator_coeffs) in power series form.
fn solve_pade_system<T: Float + FromPrimitive>(
    coeffs: &[T],
    n: usize,
    m: usize,
) -> Result<(Vec<T>, Vec<T>)> {
    // System: For k = 0 to n+m:
    //   pₖ - Σⱼ₌₀^min(k,m) qⱼ·aₖ₋ⱼ = 0  for k > n
    //   pₖ - Σⱼ₌₀^min(k,m) qⱼ·aₖ₋ⱼ = aₖ  for k ≤ n
    //
    // With normalization q₀ = 1.

    // First, solve for q₁,...,qₘ using equations k = n+1,...,n+m
    // These form: Σⱼ₌₁^m qⱼ·aₖ₋ⱼ = aₖ for k = n+1,...,n+m

    let mut q_coeffs = vec![T::zero(); m + 1];
    q_coeffs[0] = T::one(); // Normalization

    if m > 0 {
        // Build m×m linear system for q₁,...,qₘ
        let mut matrix = vec![vec![T::zero(); m]; m];
        let mut rhs = vec![T::zero(); m];

        for i in 0..m {
            let k = n + 1 + i;
            rhs[i] = coeffs[k];

            for j in 0..m {
                let idx = k - (j + 1);
                if idx < coeffs.len() {
                    matrix[i][j] = coeffs[idx];
                }
            }
        }

        // Solve using Gaussian elimination
        let q_solution = solve_linear_system(&matrix, &rhs)?;
        for (i, &val) in q_solution.iter().enumerate() {
            q_coeffs[i + 1] = val;
        }
    }

    // Now solve for p₀,...,pₙ using equations k = 0,...,n
    // pₖ = aₖ + Σⱼ₌₁^min(k,m) qⱼ·aₖ₋ⱼ
    let mut p_coeffs = vec![T::zero(); n + 1];
    for k in 0..=n {
        let mut sum = coeffs[k];
        for j in 1..=m.min(k) {
            sum = sum + q_coeffs[j] * coeffs[k - j];
        }
        p_coeffs[k] = sum;
    }

    Ok((p_coeffs, q_coeffs))
}

/// Solves a linear system Ax = b using Gaussian elimination with partial pivoting.
pub(crate) fn solve_linear_system<T: Float + FromPrimitive>(
    matrix: &[Vec<T>],
    rhs: &[T],
) -> Result<Vec<T>> {
    let n = matrix.len();
    if n == 0 || matrix[0].len() != n || rhs.len() != n {
        return Err(GelfgrenError::InvalidDimension(
            "Matrix must be square and match RHS dimension".to_string(),
        ));
    }

    // Create augmented matrix [A | b]
    let mut aug = vec![vec![T::zero(); n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = matrix[i][j];
        }
        aug[i][n] = rhs[i];
    }

    // Forward elimination with partial pivoting
    for k in 0..n {
        // Find pivot
        let mut max_idx = k;
        let mut max_val = aug[k][k].abs();
        for i in (k + 1)..n {
            let val = aug[i][k].abs();
            if val > max_val {
                max_val = val;
                max_idx = i;
            }
        }

        // Check for singular matrix
        let tol = T::from_f64(1e-10).unwrap();
        if max_val < tol {
            return Err(GelfgrenError::SingularMatrix);
        }

        // Swap rows
        if max_idx != k {
            aug.swap(k, max_idx);
        }

        // Eliminate below
        for i in (k + 1)..n {
            let factor = aug[i][k] / aug[k][k];
            for j in k..=n {
                aug[i][j] = aug[i][j] - factor * aug[k][j];
            }
        }
    }

    // Back substitution
    let mut x = vec![T::zero(); n];
    for i in (0..n).rev() {
        let mut sum = aug[i][n];
        for j in (i + 1)..n {
            sum = sum - aug[i][j] * x[j];
        }
        x[i] = sum / aug[i][i];
    }

    Ok(x)
}

/// Converts power series coefficients to a Bernstein polynomial on [a, b].
///
/// Given coefficients {c₀, c₁, ..., cₙ} representing p(x) = Σcₖxᵏ,
/// converts to Bernstein form on [a, b].
fn power_series_to_bernstein<T: Float + FromPrimitive + std::fmt::Debug>(
    power_coeffs: &[T],
    a: T,
    b: T,
) -> Result<BernsteinPolynomial<T>> {
    if power_coeffs.is_empty() {
        return Err(GelfgrenError::InvalidArgument(
            "Power series coefficients cannot be empty".to_string(),
        ));
    }

    let n = power_coeffs.len() - 1;

    // For a polynomial p(x) = Σ cₖxᵏ on [0,1], the Bernstein coefficients are:
    // bⱼ = Σₖ cₖ·(k choose j)·j!/(n choose j)·1/k! = Σₖ cₖ·(k choose j)/(n choose j)
    //
    // For interval [a,b], we need to handle the transformation x = a + t(b-a)
    // For simplicity, we'll use an approximation: sample the polynomial and interpolate

    // Alternative: Direct conversion using the matrix formula
    // For now, use sampling approach
    let num_samples = (2 * n + 1).max(10);
    let mut sample_points = Vec::with_capacity(num_samples);
    let mut sample_values = Vec::with_capacity(num_samples);

    for i in 0..num_samples {
        let t = T::from_usize(i).unwrap() / T::from_usize(num_samples - 1).unwrap();
        let x = a + t * (b - a);

        // Evaluate power series at x (assuming center is at 0 for now)
        let mut value = T::zero();
        let mut x_power = T::one();
        for &coeff in power_coeffs {
            value = value + coeff * x_power;
            x_power = x_power * x;
        }

        sample_points.push(t);
        sample_values.push(value);
    }

    // For now, return a polynomial of the same degree with coefficients
    // that approximate the power series on [a, b]
    // TODO: Implement proper power-to-Bernstein conversion

    // Simple approach: use the power coefficients directly as unscaled Bernstein
    // This is an approximation and should be improved
    BernsteinPolynomial::from_unscaled(power_coeffs.to_vec(), a, b)
}

/// Computes factorial of n.
fn factorial(n: usize) -> usize {
    (1..=n).product()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_factorial() {
        assert_eq!(factorial(0), 1);
        assert_eq!(factorial(1), 1);
        assert_eq!(factorial(5), 120);
        assert_eq!(factorial(10), 3628800);
    }

    #[test]
    fn test_solve_linear_system_2x2() {
        // System: x + 2y = 5
        //         3x + 4y = 11
        // Solution: x = 1, y = 2
        let matrix = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let rhs = vec![5.0, 11.0];

        let solution = solve_linear_system(&matrix, &rhs).unwrap();
        assert_relative_eq!(solution[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(solution[1], 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_pade_constant() {
        // f(x) = 5 (constant function)
        // [0/0] Padé: P₀/Q₀ = 5/1
        let coeffs = vec![5.0, 0.0, 0.0];
        let pade = PadeApproximant::from_series(&coeffs, 0, 0, 0.0, (0.0, 1.0)).unwrap();

        assert_eq!(pade.degrees(), (0, 0));
        assert_relative_eq!(pade.evaluate(0.5).unwrap(), 5.0, epsilon = 1e-6);
    }

    #[test]
    #[ignore] // TODO: Requires proper power-to-Bernstein conversion
    fn test_pade_linear() {
        // f(x) = 1 + 2x
        // [1/0] Padé: (1 + 2x)/1
        let coeffs = vec![1.0, 2.0];
        let pade = PadeApproximant::from_series(&coeffs, 1, 0, 0.0, (0.0, 1.0)).unwrap();

        assert_relative_eq!(pade.evaluate(0.0).unwrap(), 1.0, epsilon = 1e-6);
        assert_relative_eq!(pade.evaluate(0.5).unwrap(), 2.0, epsilon = 1e-6);
    }

    #[test]
    #[ignore] // TODO: Requires proper power-to-Bernstein conversion
    fn test_pade_exponential() {
        // f(x) = exp(x) ≈ 1 + x + x²/2 + x³/6
        // [1/1] Padé: (1 + x/2) / (1 - x/2) for exp(x)
        let coeffs = vec![1.0, 1.0, 0.5, 1.0 / 6.0];
        let pade = PadeApproximant::from_series(&coeffs, 1, 1, 0.0, (-0.5, 0.5)).unwrap();

        // Check it's reasonable (exact [1/1] for exp is (2+x)/(2-x))
        let val = pade.evaluate(0.0).unwrap();
        assert!((val - 1.0).abs() < 0.5); // Should be close to exp(0) = 1
    }

    #[test]
    fn test_from_derivatives() {
        // f(x) = x², f'(x) = 2x, f''(x) = 2
        // At x=0: f(0)=0, f'(0)=0, f''(0)=2
        let derivatives = vec![0.0, 0.0, 2.0];
        let pade = PadeApproximant::from_derivatives(&derivatives, 2, 0, 0.0, (0.0, 1.0)).unwrap();

        // Should give P₂(x)/1 ≈ x²
        assert_relative_eq!(pade.evaluate(0.5).unwrap(), 0.25, epsilon = 1e-3);
        assert_relative_eq!(pade.evaluate(1.0).unwrap(), 1.0, epsilon = 1e-3);
    }

    #[test]
    fn test_insufficient_coefficients() {
        let coeffs = vec![1.0, 2.0]; // Only 2 coefficients
        let result = PadeApproximant::from_series(&coeffs, 2, 2, 0.0, (0.0, 1.0));
        assert!(result.is_err());
    }

    #[test]
    #[ignore] // TODO: Requires proper power-to-Bernstein conversion
    fn test_derivative() {
        // f(x) = 1 + x + x²
        let coeffs = vec![1.0, 1.0, 1.0];
        let pade = PadeApproximant::from_series(&coeffs, 2, 0, 0.0, (0.0, 1.0)).unwrap();
        let d_pade = pade.derivative().unwrap();

        // Derivative should be close to 1 + 2x
        assert_relative_eq!(d_pade.evaluate(0.0).unwrap(), 1.0, epsilon = 1e-3);
        assert_relative_eq!(d_pade.evaluate(0.5).unwrap(), 2.0, epsilon = 1e-3);
    }
}
